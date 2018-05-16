import logging
from pycipapi.cipapi_client import CipApiClient
from gel2decipher.clients.decipher_client import DecipherClient
import gel2decipher.models.gel2decipher_mappings as gel2decipher
from gel2decipher.models.decipher_models import *
from protocols.migration.migration_helpers import MigrationHelpers
from protocols.cva_1_0_0 import ReportEventEntry, ObservedVariant, VariantRepresentation, Assembly, VariantAvro, \
    VariantAnnotation, ConsequenceType, VariantCall, ReportEvent, Tier, SequenceOntologyTerm, HpoTerm, TernaryOption
from requests import HTTPError
from pyark.cva_client import CvaClient
from DataModels.GelPedigree import GelRDParticipant


class UnacceptableCase(ValueError):

    pass


class Gel2Decipher(object):

    def __init__(self, config):
        Gel2Decipher._sanity_checks(config)
        self.gel_user = config['gel_user']
        self.gel_password = config['gel_password']
        self.cipapi_url = config['cipapi_url']
        self.cipapi = CipApiClient(self.cipapi_url, user=self.gel_user, password=self.gel_password)
        self.cva_url = config['cva_url']
        self.cva = CvaClient(self.cva_url, user=self.gel_user, password=self.gel_password)
        self.report_events_client = self.cva.report_events()
        self.decipher_system_key = config['decipher_system_key']
        self.decipher_user_key = config['decipher_user_key']
        self.decipher_url = config['decipher_url']
        self.send_absent_phenotypes = config['send_absent_phenotypes']
        self.decipher = DecipherClient(self.decipher_url, self.decipher_system_key, self.decipher_user_key)

    @staticmethod
    def _sanity_checks(config):
        assert config is not None, "Empty config!"
        # TODO!

    def _get_pedigree(self, case_id, case_version):
        interpretation_request_json = self.cipapi.get_interpretation_request(case_id, case_version)
        pedigree = self.cipapi.get_pedigree(interpretation_request_json)
        return pedigree

    @staticmethod
    def _get_person_id_by_relation(persons, relation):
        person_id = None
        for person in persons:
            if person['relation'] == relation:
                person_id = person['person_id']
        return person_id

    @staticmethod
    def _get_proband_observed_variant(observed_variants, proband_id, obfustcated=True):
        """
        :type observed_variants: list
        :type proband_id: str
        :rtype: ObservedVariant
        """
        proband = None
        if obfustcated:
            proband_id = gel2decipher.hash_id(proband_id)
        for observed_variant in observed_variants:  # type: ObservedVariant
            variant_call = observed_variant.variantCall     # type: VariantCall
            ov_participant_id = variant_call.participantId
            if obfustcated:
                ov_participant_id = gel2decipher.hash_id(ov_participant_id)
            if proband_id == ov_participant_id:
                proband = observed_variant
                break
        return proband

    @staticmethod
    def _get_variant_representation_grch37(observed_variant):
        """
        :type observed_variant: ObservedVariant
        :rtype: VariantAvro
        """
        grch37_variant = None
        for variant_representation in observed_variant.variant.variants:  # type: VariantRepresentation
            if variant_representation.assembly == Assembly.GRCh37:
                grch37_variant = variant_representation.variant  # type: VariantAvro
                break
        return grch37_variant

    @staticmethod
    def _select_consequence_type(so_terms, gene_symbols, tier):
        """
        :type so_terms: list
        :type gene_symbols: list
        :type tier: Tier
        :rtype: ConsequenceType
        """
        # filter consequence types by the provided list of gene symbols
        filtered_cts = [ct for ct in so_terms if ct.geneName in gene_symbols]
        if len(filtered_cts) == 0:
            filtered_cts = so_terms

        # filter consequence types by SO terms
        so_terms = {
            Tier.TIER1: ["SO:0001893", "SO:0001574", "SO:0001575", "SO:0001587", "SO:0001589", "SO:0001578",
                         "SO:0001582"],
            Tier.TIER2: ["SO:0001889", "SO:0001821", "SO:0001822", "SO:0001583", "SO:0001630", "SO:0001626"]
        }
        filtered_cts_by_so = []
        for ct in filtered_cts:  # type: ConsequenceType
            sos = set([so.accession for so in ct.sequenceOntologyTerms])  # type: SequenceOntologyTerm
            matching_sos = sos.intersection(so_terms[tier])
            if len(matching_sos) > 0:
                filtered_cts_by_so.append(ct)

        # filter consequence types by biotypes
        biotypes = ["IG_C_gene", "IG_D_gene"," IG_J_gene", "IG_V_gene", "IG_V_gene", "protein_coding",
                    "nonsense_mediated_decay", "non_stop_decay", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"]
        filtered_cts_by_bt = [ct for ct in filtered_cts_by_so if ct.biotype in biotypes]    # type: ConsequenceType

        # filter consequence types by flags
        transcript_flags = set(["basic"])
        filtered_cts_by_flag = [
            ct for ct in filtered_cts_by_bt
            if len(transcript_flags.intersection(ct.transcriptAnnotationFlags)) > 0]  # type: ConsequenceType

        # sort by transcript so the selection is deterministic
        filtered_cts_by_flag.sort(key=lambda x: x.ensemblTranscriptId)
        consequence_type = filtered_cts_by_flag[0]  # type: ConsequenceType
        logging.info("The selected transcript is {} at gene {}".format(
            consequence_type.ensemblTranscriptId, consequence_type.geneName))

        return consequence_type

    def _send_pedigree_member_phenotypes(self, pedigree_member, decipher_person_id):
        """
        :type pedigree_member: PedigreeMember
        :type decipher_person_id: str
        :return:
        """
        accepted_phenotypes = []
        rejected_phenotypes = []
        decipher_phenotype_ids = []
        for phenotype in pedigree_member.hpoTermList:   # type: HpoTerm
            # avoid sending unknown presence phenotypes
            if phenotype.termPresence == TernaryOption.unknown:
                logging.warn("Skipping phenotype {}".format(phenotype.term))
                continue
            # optionally send absent phenotypes
            elif not self.send_absent_phenotypes and phenotype.termPresence == TernaryOption.no:
                logging.warn("Skipping phenotype {}".format(phenotype.term))
                continue
            else:
                dec_phenotype = gel2decipher.map_phenotype(phenotype, decipher_person_id)
                try:
                    phenotype_id = self.decipher.create_phenotypes([dec_phenotype], decipher_person_id)[0]
                    decipher_phenotype_ids.append(phenotype_id)
                    accepted_phenotypes.append(phenotype)
                except HTTPError:
                    logging.warning("Rejected phenotype: {}".format(phenotype.toJsonString()))
                    rejected_phenotypes.append(phenotype)
        return accepted_phenotypes, rejected_phenotypes, decipher_phenotype_ids

    def send_case(self, case_id, case_version):

        # fetch pedigree and get proband
        case = self.cipapi.get_case(case_id, case_version)
        pedigree = case.get_pedigree()
        proband = pedigree.get_proband()
        father = pedigree.get_father(proband)   # type: GelRDParticipant
        mother = pedigree.get_mother(proband)   # type: GelRDParticipant
        logging.info("The proband is {}".format(proband.participantId))

        # fetch variants from CVA
        report_events_iterator = self.report_events_client.get_report_events({
            'parent_id': case_id, 'parent_version': case_version, 're_type': 'tiered', 'tier': 'TIER1,TIER2',
            'vcf_format': True
        })

        # filters down the variants
        accepted_variants = []
        for report_event in report_events_iterator:  # type: ReportEventEntry
            # selects the observed variant for the proband
            proband_ov = Gel2Decipher._get_proband_observed_variant(
                report_event.observedVariants, proband_id=proband.participantId)  # type: ObservedVariant
            if proband_ov is None:
                raise UnacceptableCase("There is a report event with no observed variant for the proband")
            variant_call = proband_ov.variantCall  # type: VariantCall

            # selects the variant representation for assembly GRCh37
            grch37_variant = Gel2Decipher._get_variant_representation_grch37(proband_ov)
            if grch37_variant is None:
                logging.warning("The report event does not have coordinates in GRCh37")
                continue
            else:
                accepted_variants.append((grch37_variant, report_event, variant_call))

        if len(accepted_variants) == 0:
            message = "The case id={} and version={} has no variants".format(case_id, case_version)
            logging.warning(message)
            raise UnacceptableCase(message)

        # create patient for proband in Decipher
        try:
            patient_id = self.decipher.create_patients(
                [gel2decipher.map_pedigree_member_to_patient(
                    proband, self.decipher.project_id, self.decipher.user_id)])[0]['patient_id']
        except HTTPError:
            # the patient already exists??
            raise UnacceptableCase("Patient registration failed, don't know how to continue")

        # create phenotypes
        dec_persons = self.decipher.get_persons_by_patient(patient_id)
        logging.info("The persons: " + str(dec_persons))
        dec_proband = Gel2Decipher._get_person_id_by_relation(dec_persons, 'patient')
        dec_mother = Gel2Decipher._get_person_id_by_relation(dec_persons, 'mother')
        dec_father = Gel2Decipher._get_person_id_by_relation(dec_persons, 'father')

        self._send_pedigree_member_phenotypes(proband, dec_proband)

        # updates mother and father affection status because they are created automatically
        if mother is not None:
            self.decipher.update_person(gel2decipher.map_affection_status(mother.affectionStatus), dec_mother)
            self._send_pedigree_member_phenotypes(mother, dec_mother)
        if father is not None:
            self.decipher.update_person(gel2decipher.map_affection_status(father.affectionStatus), dec_father)
            self._send_pedigree_member_phenotypes(father, dec_father)

        for member in pedigree.members:   # type: GelRDParticipant
            logging.info("Member of family: {}".format(member.pedigreeId))
            if member.pedigreeId not in [proband.pedigreeId, father.pedigreeId, mother.pedigreeId]:
                dec_person = self.decipher.create_persons(
                    [gel2decipher.map_pedigree_member_to_person(
                        member, patient_id, pedigree.get_relationship(member.pedigreeId, proband.pedigreeId))],
                    patient_id)[0]
                self._send_pedigree_member_phenotypes(member, dec_person)

        # push the variants to Decipher
        unique_variants = {}
        for variant, report_event, variant_call in accepted_variants:  # type: (VariantAvro, ReportEventEntry)
            annotation = variant.annotation  # type: VariantAnnotation
            gene_symbols = [x.geneSymbol for x in report_event.reportEvent.genomicEntities]
            consequence_type = Gel2Decipher._select_consequence_type(
                annotation.consequenceTypes, gene_symbols, report_event.reportEvent.tier)  # type: ConsequenceType

            # builds the variant in decipher model
            dec_variant = gel2decipher.map_report_event(
                report_event.reportEvent, variant, variant_call, consequence_type, patient_id)
            uid = "{}:{}:{}:{}".format(dec_variant.chr, dec_variant.start, dec_variant.ref_allele,
                                       dec_variant.alt_allele)
            if uid not in unique_variants:
                unique_variants[uid] = dec_variant
        # NOTE: this removes the duplicated variants from composite heterozygous report events
        self.decipher.create_snvs(unique_variants.values(), patient_id)

        return patient_id

import logging
from gel2decipher.models.decipher_models import *
from protocols.participant_1_0_3 import PedigreeMember, PersonKaryotipicSex
from protocols.cva_1_0_0 import VariantAvro, VariantCall, ConsequenceType, Assembly
import hashlib
from datetime import datetime


def map_sex(gel_sex):
    sex_map = {
        "female": "46XX",
        "male": "46XY",
        "unknown": "unknown",
        "undetermined": "unknown"
    }
    return sex_map.get(gel_sex, None)


def map_kariotypic_sex(gel_kariotypic_sex):
    """
    When mapping to Decipher we lose XXYY, XXXY and XXXX.
    :param gel_kariotypic_sex: 
    :return: 
    """
    keriotypic_sex_map = {
        "XX": "46XX",
        "XY": "46XY",
        "XO": "45X",
        "XXY": "47XXY",
        "XXX": "47XX",
        "XYY": "47XYY",
        "XXYY": "other",
        "XXXY": "other",
        "XXXX": "other",
        "OTHER": "other",
        "UNKNOWN": "unknown"
    }
    return keriotypic_sex_map.get(gel_kariotypic_sex, "unknown")


def map_patient(participant, project_id, user_id):
    """

    :param user_id:
    :param project_id:
    :param participant:
    :type participant: RDParticipant
    :return:
    """
    # TODO: parse kariotypic_sex = proband.personKaryotipicSex
    patient = Patient(
        sex=map_sex(participant.sex),
        reference=hash_id(participant.gelId),
        project_id=project_id,
        consent="No",
        user_id=user_id
    )
    logging.warn("Creating patient with id={}".format(patient.reference))
    return patient


def map_yob_to_age(yob):
    """
    We assume it is the age when the patient was registered
    :param yob:
    :return:
    """
    if yob is not None and yob > 0:
        today = datetime.today()
        age = str(today.year - yob)
    else:
        age = 'unknown'
    return age


def obfuscate_pedigree_member(member):
    member.participantId = hash_id(member.participantId)
    member.pedigreeId = hash_id(member.pedigreeId)
    member.gelSuperFamilyId = hash_id(member.gelSuperFamilyId)
    member.fatherId = hash_id(member.fatherId)
    member.motherId = hash_id(member.motherId)
    member.superFatherId = hash_id(member.superFatherId)
    member.superMotherId = hash_id(member.superMotherId)
    return member


def map_pedigree_member_to_patient(member, project_id, user_id):
    """

    :param user_id:
    :param project_id:
    :param member:
    :type member: PedigreeMember
    :return:
    """

    member = obfuscate_pedigree_member(member)
    patient = Patient(
        sex=map_kariotypic_sex(member.personKaryotypicSex),
        reference=member.participantId,
        project_id=project_id,
        user_id=user_id,
        age=map_yob_to_age(member.yearOfBirth),
        # we mark aneuploidy only for sexual chromosomes, otherwise None as we don't have evidence for other chromosomes
        aneuploidy=True if member.personKaryotypicSex in [
            PersonKaryotipicSex.XO, PersonKaryotipicSex.XXX, PersonKaryotipicSex.XXXX, PersonKaryotipicSex.XXXY,
            PersonKaryotipicSex.XXY, PersonKaryotipicSex.XXYY, PersonKaryotipicSex.XYY
        ] else None,
        # we are missing the carrier status consent field
        consent='Yes' if member.consentStatus.secondaryFindingConsent else 'No',
        note="\n".join(["{term}({presence})".format(term=hpo.term, presence=hpo.termPresence)
                        for hpo in member.hpoTermList])
    )
    return patient


def map_phenotype(phenotype, person_id):
    map_presence = {
        "yes": "present",
        "no": "absent"
    }
    decipher_phenotype = Phenotype(
        person_id=person_id,
        phenotype_id=int(phenotype.term.replace("HP:", ""))
    )
    observation = map_presence.get(phenotype.termPresence, None)
    if observation:
        decipher_phenotype.observation = observation
    logging.info("Creating phenotype {}:{}".format(
        decipher_phenotype.phenotype_id, decipher_phenotype.observation))
    return decipher_phenotype


def map_variant(variant, patient_id, gel_id):
    """

    :param gel_id:
    :param variant:
    :type variant: protocols.reports_5_0_0.ReportedVariant
    :param patient_id:
    :return:
    """
    variant_call = None
    for vc in variant.variantCalls:
        if vc.participantId == gel_id:
            variant_call = vc
    if variant_call is None:
        raise ValueError("No called genotype for the proband")
    report_event = variant.reportEvents[0]
    gene_symbol = report_event.genomicEntities[0].geneSymbol
    snv = Snv(
        patient_id=patient_id,
        assembly=map_assembly(variant.variantCoordinates.assembly),
        chr=normalise_chromosome(variant.variantCoordinates.chromosome),
        start=variant.variantCoordinates.position,
        ref_allele=variant.variantCoordinates.reference,
        alt_allele=variant.variantCoordinates.alternate,
        genotype=map_genotype(variant_call.zygosity),
        intergenic=True,    # this is a temp solution while we figure out how to fetch a correct transcript
        #user_transcript=None,
        user_gene=gene_symbol
    )
    # TODO: map user transcript
    # TODO: map user gene
    # TODO: map inheritance from report event

    return snv


def map_report_event(grch37_variant, variant_call, consequence_type, patient_id):
    """

    :param grch37_variant:
    :type grch37_variant: VariantAvro
    :param variant_call:
    :type variant_call: VariantCall
    :param consequence_type:
    :type consequence_type: ConsequenceType
    :return:
    """
    snv = Snv(
        patient_id=patient_id,
        assembly=map_assembly(Assembly.GRCh37),
        chr=normalise_chromosome(grch37_variant.chromosome),
        start=grch37_variant.start,
        ref_allele=grch37_variant.reference,
        alt_allele=grch37_variant.alternate,
        genotype=map_genotype(variant_call.zygosity),
        intergenic=False,
        user_transcript=consequence_type.ensemblTranscriptId,
        user_gene=consequence_type.geneName
    )
    logging.warn("Creating variant {}:{}:{}:{}:{}".format(
        snv.assembly, snv.chr, snv.start, snv.ref_allele, snv.alt_allele))

    return snv


def map_genotype(gel_genotype):

    genotype_map = {
        "reference_homozygous": None,
        "heterozygous": "Heterozygous",
        "alternate_homozygous": "Homozygous",
        "missing": None,
        "half_missing_reference": None,
        "half_missing_alternate": None,
        "alternate_hemizygous": "Hemizygous",
        "reference_hemizygous": None,
        "unk": None
    }
    return genotype_map.get(gel_genotype, None)


def map_assembly(gel_assembly):
    assembly_map = {
        "GRCh37": "GRCh37/hg19",
        "GRCh38": None  # NOTE: this is not supported
    }
    return assembly_map.get(gel_assembly, None)


def normalise_chromosome(chromosome):
    return chromosome.replace("chr", "")


def hash_id(identifier):
    return hashlib.sha224(str(identifier)).hexdigest()

from gel2decipher.models.decipher_models import *
from protocols.reports_3_0_0 import RDParticipant
import hashlib


def map_sex(gel_sex):
    sex_map = {
        "female": "46XX",
        "male": "46XY",
        "unknown": "unknown",
        "undetermined": "unknown"
    }
    return sex_map.get(gel_sex, None)


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
    return patient


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
    return hashlib.sha224(identifier).hexdigest()

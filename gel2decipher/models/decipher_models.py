from booby import Model, fields
from enum import Enum


class Patient(Model):
    # required fields
    sex = fields.Field(choices=['46XY', '46XX', '45X', '47XX', '47XXY', '47XYY', 'other', 'unknown'], required=True)
    reference = fields.String(required=True)
    project_id = fields.Integer(required=True)
    # optional fields
    age = fields.Field(choices=['unknown', 'Prenatal'] + [str(x) for x in range(0, 101)], default="unknown")
    prenatal = fields.Field(choices=[x for x in range(10, 43)].append(None))
    aneuploidy = fields.Boolean(default=False)
    user_id = fields.Integer()
    note = fields.String()
    consent = fields.Field(choices=['No', 'Yes'], default="No")


class Assembly(Enum):
    grch37 = "GRCh37/hg19"


class Genotype(Enum):
    homozygous = "Homozygous"
    heterozygous = "Heterozygous"
    hemizygous = "Hemizygous"
    mitochondrial_homoplasmy = "Mitochondrial Homoplasmy",
    mitochondrial_heteroplasmy = "Mitochondrial Heteroplasmy"


class Inheritance(Enum):
    unknown = "Unknown"
    de_novo_constitutive = "De novo constitutive"
    de_novo_mosaic = "De novo mosaic",
    paternally_constitutive = "Paternally inherited, constitutive in father"
    paternally_mosaic = "Paternally inherited, mosaic in father"
    maternally_constitutive = "Maternally inherited, constitutive in mother"
    maternally_mosaic = "Maternally inherited, mosaic in mother"
    biparental = "Biparental"
    imbalance = "Imbalance arising from a balanced parental rearrangement"


class ClinicalSignificance(Enum):
    pathogenic = "Pathogenic"
    likely_pathogenic = "Likely pathogenic"
    uncertain = "Uncertain"
    likely_benign = "Likely benign"
    benign = "Benign"


class Penetrance(Enum):
    full = "Full"
    partial = "Partial"
    uncertain = "Uncertain"
    none = "None"


class Snv(Model):
    # required fields
    patient_id = fields.Integer(required=True)
    assembly = fields.Field(choices=[x.value for x in Assembly.__members__.values()], required=True)
    chr = fields.Field(choices=["X", "Y", "MT"] + [str(x) for x in range(1, 23)], required=True)
    start = fields.Integer(required=True)
    ref_allele = fields.String(required=True)
    alt_allele = fields.String(required=True)
    genotype = fields.Field(choices=[x.value for x in Genotype.__members__.values()], required=True)
    # optional fields
    user_transcript = fields.String()
    user_gene = fields.String()
    intergenic = fields.Boolean(default=False)
    inheritance = fields.Field(choices=[x.value for x in Inheritance.__members__.values()],
                               default=Inheritance.unknown.value)
    pathogenicity = fields.Field(choices=[x.value for x in ClinicalSignificance.__members__.values()].append(None))
    contribution = fields.Field(choices=[x.value for x in Penetrance.__members__.values()].append(None))
    shared = fields.String()


class Cnv(Model):
    # TODO
    pass


class Relation(Enum):
    father = "father"
    mother = "mother"
    brother = "brother"
    sister = "sister"
    son = "son"
    daughter = "daughter"
    maternal_half_brother = "maternal_half_brother"
    maternal_half_sister = "maternal_half_sister"
    maternal_grandmother = "maternal_grandmother"
    maternal_grandfather = "maternal_grandfather"
    maternal_aunt = "maternal_aunt"
    maternal_uncle = "maternal_uncle"
    maternal_cousin = "maternal_cousin"
    paternal_half_brother = "paternal_half_brother"
    paternal_half_sister = "paternal_half_sister"
    paternal_grandmother = "paternal_grandmother"
    paternal_grandfather = "paternal_grandfather"
    paternal_aunt = "paternal_aunt"
    paternal_uncle = "paternal_uncle"
    paternal_cousin = "paternal_cousin"
    fraternal_nephew = "fraternal_nephew"
    fraternal_niece = "fraternal_niece"
    sororal_nephew = "sororal_nephew"
    sororal_niece = "sororal_niece"
    grandson_through_daughter = "grandson_through_daughter"
    grandson_through_son = "grandson_through_son"
    granddaughter_through_daughter = "granddaughter_through_daughter"
    granddaughter_through_son = "granddaughter_through_son"
    other_maternal_relative = "other_maternal_relative"
    other_paternal_relative = "other_paternal_relative"
    other_blood_relative = "other_blood_relative"


class AffectionStatus(Enum):
    affected = "affected"
    unaffected = "unaffected"
    unknown = "unknown"


class Person(Model):
    # required fields
    patient_id = fields.Integer(required=True)
    relation = fields.Field(
        choices=filter(lambda x: x not in [Relation.father.value, Relation.mother.value],
                       [x.value for x in Relation.__members__.values()]),
        required=True)
    # optional fields
    relation_status = fields.Field(choices=[x.value for x in AffectionStatus.__members__.values()],
                                   default=AffectionStatus.unknown.value)


class TermPresence(Enum):
    present = "present"
    absent = "absent"


class Phenotype(Model):
    # required fields
    person_id = fields.Integer(required=True)
    phenotype_id = fields.Integer(required=True)
    # optional fields
    observation = fields.Field(choices=[x.value for x in TermPresence.__members__.values()],
                               default=TermPresence.present.value)


class SnvPathogenicity(Model):
    # TODO
    pass

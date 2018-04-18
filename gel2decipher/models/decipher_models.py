from booby import Model, fields


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


class Snv(Model):
    # required fields
    patient_id = fields.Integer(required=True)
    assembly = fields.Field(choices=["GRCh37/hg19"], required=True)
    chr = fields.Field(choices=["X", "Y", "MT"] + [str(x) for x in range(1, 23)], required=True)
    start = fields.Integer(required=True)
    ref_allele = fields.String(required=True)
    alt_allele = fields.String(required=True)
    genotype = fields.Field(choices=["Homozygous", "Heterozygous", "Hemizygous", "Mitochondrial Homoplasmy",
                                     "Mitochondrial Heteroplasmy"], required=True)
    # optional fields
    user_transcript = fields.String()
    user_gene = fields.String()
    intergenic = fields.Boolean(default=False)
    inheritance = fields.Field(choices=["Unknown", "De novo constitutive", "De novo mosaic",
                                        "Paternally inherited, constitutive in father",
                                        "Paternally inherited, mosaic in father",
                                        "Maternally inherited, constitutive in mother",
                                        "Maternally inherited, mosaic in mother", "Biparental",
                                        "Imbalance arising from a balanced parental rearrangement"],
                               default="Unknown")
    pathogenicity = fields.Field(choices=["Definitely pathogenic", "Likely pathogenic", "Uncertain", "Likely benign",
                                          "Benign", None])
    contribution = fields.Field(choices=["Full", "Partial", "Uncertain", "None", None])
    shared = fields.String()


class Cnv(Model):
    # TODO
    pass


class Person(Model):
    # required fields
    patient_id = fields.Integer(required=True)
    relation = fields.Field(choices=["brother", "sister", "son", "daughter", "maternal_half_brother",
                                     "maternal_half_sister", "maternal_grandmother", "maternal_grandfather",
                                     "maternal_aunt", "maternal_uncle", "maternal_cousin", "paternal_half_brother",
                                     "paternal_half_sister", "paternal_grandmother", "paternal_grandfather",
                                     "paternal_aunt", "paternal_uncle", "paternal_cousin", "fraternal_nephew",
                                     "fraternal_niece", "sororal_nephew", "sororal_niece", "grandson_through_daughter",
                                     "grandson_through_son", "granddaughter_through_daughter",
                                     "granddaughter_through_son", "other_maternal_relative", "other_paternal_relative",
                                     "other_blood_relative"], required=True)
    # optional fields
    relation_status = fields.Field(choices=["affected", "unaffected", "unknown"], default="unknown")


class Phenotype(Model):
    # required fields
    person_id = fields.Integer(required=True)
    phenotype_id = fields.Integer(required=True)
    # optional fields
    observation = fields.Field(choices=["present", "absent"], default="present")


class SnvPathogenicity(Model):
    # TODO
    pass

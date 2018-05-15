import os
import logging
from unittest import TestCase
from requests.exceptions import HTTPError, InvalidSchema

from gel2decipher.clients.decipher_client import DecipherClient
from gel2decipher.models.decipher_models import *
from gel2decipher.case_sender import Gel2Decipher, UnacceptableCase


class TestGel2Decipher(TestCase):
    # credentials
    SYSTEM_KEY = os.getenv("DECIPHER_SYSTEM_KEY")
    USER_KEY = os.getenv("DECIPHER_USER_KEY")
    AUTH_HEADERS = {'X-Auth-Token-System': SYSTEM_KEY, 'X-Auth-Token-User': USER_KEY}
    DECIPHER_URL_BASE = "https://decipher.sanger.ac.uk/API/"

    CIPAPI_URL_BASE = os.getenv("CIPAPI_URL")
    CVA_URL_BASE = os.getenv("CVA_URL")
    GEL_USER = os.getenv("GEL_USER")
    GEL_PASSWORD = os.getenv("GEL_PASSWORD")

    def setUp(self):
        logging.basicConfig(level=logging.INFO)
        config = {
            'gel_user': self.GEL_USER,
            'gel_password': self.GEL_PASSWORD,
            'cipapi_url': self.CIPAPI_URL_BASE,
            'cva_url': self.CVA_URL_BASE,
            'decipher_system_key': self.SYSTEM_KEY,
            'decipher_user_key': self.USER_KEY,
            'decipher_url': self.DECIPHER_URL_BASE,
            'send_absent_phenotypes': False
        }
        self.sender = Gel2Decipher(config)

    def test_send_case(self):

        cases = [("615", "1"), ("502", "1"), ("1026", "1"), ("2146", "1"), ("11148", "1"), ("9007", "1"),
                 ("333", "1"), ("3507", "1"), ("7139", "1"), ("11365", "1"), ]
        for case_id, case_version in cases:
            logging.warn("Trying to send case {} version {}".format(case_id, case_version))
            try:
                patient_id = self.sender.send_case(case_id, case_version)
                variants = self.sender.decipher.get_snvs(patient_id)
                #logging.warn("Rejected phenotypes: {}".format(str(rejected_phenotypes)))
                logging.warn(variants['snvs'])
            except UnacceptableCase, e:
                logging.error(e.message)


class TestDecipherApi(TestCase):

    # credentials
    SYSTEM_KEY = os.getenv("DECIPHER_SYSTEM_KEY")
    USER_KEY = os.getenv("DECIPHER_USER_KEY")
    AUTH_HEADERS = {'X-Auth-Token-System': SYSTEM_KEY, 'X-Auth-Token-User': USER_KEY}
    URL_BASE = "https://decipher.sanger.ac.uk/API/"

    def setUp(self):

        self.decipher = DecipherClient(url_base=self.URL_BASE, system_key=self.SYSTEM_KEY, user_key=self.USER_KEY)

        # PRE: a project and adequate users have been created
        # Step 1: create a patient

        self.patient1 = Patient(
            sex="46XY",
            reference="evaristo",
            project_id=self.decipher.project_id,
            consent="No",   # TODO: consent is set to yes, clarify if we want to share variants or not
            user_id=self.decipher.user_id,
            #  TODO: the user id must be set for consent=yes and it must be clinician/coordinator in this project
            note="this is just testing...",
            age="Prenatal",
            prenatal=10
        )
        self.patient2 = Patient(
            sex="46XX",
            reference="hermenegildo",
            project_id=self.decipher.project_id,
            consent="No",  # TODO: consent is set to yes, clarify if we want to share variants or not
            user_id=self.decipher.user_id,
            #  TODO: the user id must be set for consent=yes and it must be clinician/coordinator in this project
            note="this is just testing...",
            age="unknown"
        )

        response = self.decipher.create_patients([self.patient1, self.patient2])
        print response
        patient_ids = [x['patient_id'] for x in response]
        self.patient1_id = patient_ids[0]
        self.patient2_id = patient_ids[1]

    def tearDown(self):
        self.decipher.delete_patient(self.patient1_id)
        self.decipher.delete_patient(self.patient2_id)

    def test_existing_patient(self):
        # ... if the patient already exists
        try:
            self.decipher.create_patients([self.patient1])
            self.assertTrue(False)
        except HTTPError, ex:
            self.assertEqual(ex.response.status_code, 400, "Expected 400 for existing patient")

    def test_create_snv(self):
        # rs397507564
        variant = Snv(
            patient_id=self.patient1_id,
            assembly="GRCh37/hg19",
            chr="7",
            start=117119258,
            ref_allele="TCTC",
            alt_allele="T",
            genotype="Homozygous",
            user_transcript="ENST00000546407",
            user_gene="CFTR",
            shared="no"
        )
        variant_ids = self.decipher.create_snvs([variant], self.patient1_id)
        self.decipher.delete_snv(variant_ids[0])

    def test_duplicated_snv(self):
        # rs397507564
        variant = Snv(
            patient_id=self.patient1_id,
            assembly="GRCh37/hg19",
            chr="7",
            start=117119258,
            ref_allele="TCTC",
            alt_allele="T",
            genotype="Homozygous",
            user_transcript="ENST00000546407",
            user_gene="CFTR",
            shared="no"
        )
        variant2 = variant
        # tries to create the same variant twice in one request
        try:
            self.decipher.create_snvs([variant, variant2], self.patient1_id)
            self.assertTrue(False)
        except HTTPError, ex:
            self.assertEqual(ex.response.status_code, 400, "Expected 400 for duplicate records")

        # tries to create the same variant twice in two request
        variant1_id = self.decipher.create_snvs([variant], self.patient1_id)[0]
        try:
            self.decipher.create_snvs([variant2], self.patient1_id)[0]
            self.assertTrue(False)
        except HTTPError, ex:
            self.assertEqual(ex.response.status_code, 400, "Expected 400 for duplicate records")

        self.decipher.delete_snv(variant1_id)

    def test_variant_left_align(self):
        # rs864622168
        variant2 = Snv(
            patient_id=self.patient1_id,
            assembly="GRCh37/hg19",
            chr="1",
            start=145415387,  # this is left aligned to 145415368
            ref_allele="G",
            alt_allele="GAGG",
            genotype="Heterozygous",
            user_transcript="ENST00000336751",
            user_gene="HFE2",
            shared="no"
        )
        variant_ids = self.decipher.create_snvs([variant2], self.patient1_id)
        self.decipher.delete_snv(variant_ids[0])

    def test_wrong_chromosome(self):
        variant = Snv(
            patient_id=self.patient1_id,
            assembly="GRCh37/hg19",
            chr="chr7",
            start=117119258,
            ref_allele="TCTC",
            alt_allele="T",
            genotype="Homozygous",
            user_transcript="ENST00000546407",
            user_gene="CFTR",
            shared="no"
        )
        try:
            self.decipher.create_snvs([variant], self.patient1_id)
            self.assertTrue(False)
        except InvalidSchema, ex:
            self.assertTrue(True, "Expected invalid schema")

    def test_not_matching_reference(self):
        # deletion TGAG > T, but at that positions it should be TCTC > T
        variant = Snv(
            patient_id=self.patient1_id,
            assembly="GRCh37/hg19",
            chr="7",
            start=117119258,
            ref_allele="TGAG",
            alt_allele="T",
            genotype="Homozygous",
            user_transcript="ENST00000546407",
            user_gene="CFTR",
            shared="no"
        )
        try:
            variant_ids = self.decipher.create_snvs([variant], self.patient1_id)
            self.decipher.delete_snv(variant_ids[0])
            self.assertTrue(False)
        except HTTPError, ex:
            self.assertEqual(ex.response.status_code, 400, "Expected 400 for mismatching reference")

    def test_wrong_transcript(self):
        # the fake transcript ENST123456 does not overlap this variant
        variant = Snv(
            patient_id=self.patient1_id,
            assembly="GRCh37/hg19",
            chr="7",
            start=117119258,
            ref_allele="TCTC",
            alt_allele="T",
            genotype="Homozygous",
            user_transcript="ENST123456",
            user_gene="CFTR",
            shared="no"
        )
        try:
            variant_ids = self.decipher.create_snvs([variant], self.patient1_id)
            self.decipher.delete_snv(variant_ids[0])
            self.assertTrue(False)
        except HTTPError, ex:
            self.assertEqual(ex.response.status_code, 400, "Expected 400 for wrong transcript")

    def test_wrong_gene(self):
        # the fake transcript ENST123456 does not overlap this variant
        variant = Snv(
            patient_id=self.patient1_id,
            assembly="GRCh37/hg19",
            chr="7",
            start=117119258,
            ref_allele="TCTC",
            alt_allele="T",
            genotype="Homozygous",
            user_transcript="ENST00000546407",
            user_gene="OTHER",
            shared="no"
        )
        try:
            variant_ids = self.decipher.create_snvs([variant], self.patient1_id)
            self.decipher.delete_snv(variant_ids[0])
            self.assertTrue(False)
        except HTTPError, ex:
            self.assertEqual(ex.response.status_code, 400, "Expected 400 for wrong gene")

    def test_coordinate_overflow(self):
        # the fake transcript ENST123456 does not overlap this variant
        variant = Snv(
            patient_id=self.patient1_id,
            assembly="GRCh37/hg19",
            chr="7",
            start=1171192580000000,
            ref_allele="TCTC",
            alt_allele="T",
            genotype="Homozygous",
            user_transcript="ENST00000546407",
            user_gene="CFTR",
            shared="no"
        )
        try:
            variant_ids = self.decipher.create_snvs([variant], self.patient1_id)
            self.decipher.delete_snv(variant_ids[0])
            self.assertTrue(False)
        except HTTPError, ex:
            self.assertEqual(ex.response.status_code, 400, "Expected 400 for unexisting position")

    def test_not_existing_patient(self):
        # the fake transcript ENST123456 does not overlap this variant
        variant = Snv(
            patient_id=-12345,
            assembly="GRCh37/hg19",
            chr="7",
            start=117119258,
            ref_allele="TCTC",
            alt_allele="T",
            genotype="Homozygous",
            user_transcript="ENST00000546407",
            user_gene="CFTR",
            shared="no"
        )
        try:
            variant_ids = self.decipher.create_snvs([variant], -12345)
            self.decipher.delete_snv(variant_ids[0])
            self.assertTrue(False)
        except HTTPError, ex:
            self.assertEqual(ex.response.status_code, 404, "Expected 404 for unexisting patient")

    def test_batch_variants(self):
        pass

    def test_create_person(self):
        """
        Creates a person related to a patient
        :return:
        """
        person = Person(
            patient_id=self.patient1_id,
            relation="brother",
            relation_status="unaffected"
        )
        person2 = Person(
            patient_id=self.patient1_id,
            relation="sister",
            relation_status="unaffected"
        )
        person_ids = self.decipher.create_persons([person, person2], self.patient1_id)

        for person_id in person_ids:
            self.decipher.delete_person(person_id)

    def test_create_person_different_patients(self):
        """
        Creates a person related to a patient
        :return:
        """
        person = Person(
            patient_id=self.patient1_id,
            relation="brother",
            relation_status="unaffected"
        )
        try:
            self.decipher.create_persons([person], self.patient2_id)
            self.assertTrue(False)
        except HTTPError, ex:
            self.assertEqual(ex.response.status_code, 400, "Expected 400 for wrong patient id")

    def test_create_phenotypes(self):
        """
        Creates a phenotype related to a patient
        :return:
        """
        person = Person(
            patient_id=self.patient1_id,
            relation="brother",
            relation_status="unaffected"
        )
        person_ids = self.decipher.create_persons([person], self.patient1_id)
        phenotype = Phenotype(
            person_id=person_ids[0],
            phenotype_id=4322,
            observation="present"
        )
        phenotype_ids = self.decipher.create_phenotypes([phenotype], person_ids[0])
        self.decipher.delete_phenotype(phenotype_ids[0])
        # assigns phenotype to a patient
        persons = self.decipher.get_persons_by_patient(self.patient1_id)
        proband = None
        for person in persons:
            if person['relation'] == 'patient':
                proband = person
        phenotype = Phenotype(
            person_id=proband['person_id'],
            phenotype_id=4322,
            observation="present"
        )
        phenotype_ids = self.decipher.create_phenotypes([phenotype], proband['person_id'])
        self.decipher.delete_phenotype(phenotype_ids[0])
        self.decipher.delete_person(person_ids[0])

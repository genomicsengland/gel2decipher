import requests
import json
import os
from unittest import TestCase


class TestDecipherApi(TestCase):

    # credentials
    SYSTEM_KEY = os.getenv("DECIPHER_SYSTEM_KEY")
    USER_KEY = os.getenv("DECIPHER_USER_KEY")
    AUTH_HEADERS = {'X-Auth-Token-System': SYSTEM_KEY, 'X-Auth-Token-User': USER_KEY}
    URL_BASE = "https://decipher.sanger.ac.uk/API/"

    def setUp(self):
        # queries for the basic information about user and project
        response = requests.get("{url_base}info".format(url_base=self.URL_BASE), headers=self.AUTH_HEADERS)
        self.assertTrue(response.status_code == 200, response.content)
        self.project_id = response.json()["user"]["project"]["project_id"]
        self.user_id = response.json()["user"]["user_id"]

        # PRE: a project and adequate users have been created

        # Step 1: create a patient
        self.patient1 = {
            "sex": "46XY",
            "reference": "martin!",
            "project_id": self.project_id,
            "consent": "no",   # TODO: consent is set to yes, clarify if we want to share variants or not
            # "user_id": user_id,
            #  TODO: the user id must be set for consent=yes and it must be clinician/coordinator in this project
            "note": "this is just testing...",
            "age": "Prenatal",
            "prenatal": 10
        }
        self.patient2 = {
            "sex": "46XX",
            "reference": "claudia!",
            "project_id": self.project_id,
            "consent": "no",  # TODO: consent is set to yes, clarify if we want to share variants or not
            # "user_id": user_id,
            #  TODO: the user id must be set for consent=yes and it must be clinician/coordinator in this project
            "note": "this is just testing...",
            "age": "unknown"
        }
        response = requests.post("{url_base}projects/{project_id}/patients".format(
            url_base=self.URL_BASE, project_id=self.project_id),
            headers=self.AUTH_HEADERS, data=json.dumps([self.patient1, self.patient2]))
        self.assertTrue(response.status_code == 200, response.content)
        patient_ids = [x["patient_id"] for x in response.json()]
        self.patient1_id = patient_ids[0]
        self.patient2_id = patient_ids[1]

    def tearDown(self):
        response = requests.delete("{url_base}patients/{patient_id}".format(
            url_base=self.URL_BASE, patient_id=self.patient1_id), headers=self.AUTH_HEADERS)
        self.assertTrue(response.status_code == 200, response.content)
        response = requests.delete("{url_base}patients/{patient_id}".format(
            url_base=self.URL_BASE, patient_id=self.patient2_id), headers=self.AUTH_HEADERS)
        self.assertTrue(response.status_code == 200, response.content)

    def test_existing_patient(self):
        # ... if the patient already exists
        response = requests.post("{url_base}projects/{project_id}/patients".format(
            url_base=self.URL_BASE, project_id=self.project_id),
            headers=self.AUTH_HEADERS, data=json.dumps([self.patient1]))
        self.assertTrue(response.status_code == 400)
        self.assertTrue(response.json()["errors"][0] ==
                        "A patient with this internal reference number or id already exists for this project.")

    def test_variant1(self):
        # rs397507564
        variant = {
            "patient_id": self.patient1_id,
            "assembly": "GRCh37/hg19",
            "chr": "7",
            "start": 117119258,
            "ref_allele": "TCTC",
            "alt_allele": "T",
            "genotype": "Homozygous",
            "user_transcript": "ENST00000546407",
            "user_gene": "CFTR",
            "shared": "no"
        }
        variants = self.create_variants([variant], self.patient1_id)
        variant_ids = [x["patient_snv_id"] for x in variants]
        for x in variants:
            self.assertTrue(x["patient_snv_id"] is not None, "Empty snv id?")
        self.delete_variants(variant_ids)

    def test_variant_left_align(self):
        # rs864622168
        variant2 = {
            "patient_id": self.patient1_id,
            "assembly": "GRCh37/hg19",
            "chr": "1",
            "start": 145415387,  # this is left aligned to 145415368
            "ref_allele": "G",
            "alt_allele": "GAGG",
            "genotype": "Heterozygous",
            "user_transcript": "ENST00000336751",
            "user_gene": "HFE2",
            "shared": "no"
        }
        variants = self.create_variants([variant2], self.patient1_id)
        variant_ids = [x["patient_snv_id"] for x in variants]
        for x in variants:
            self.assertTrue(x["patient_snv_id"] is not None, "Empty snv id?")
        self.delete_variants(variant_ids)

    def test_wrong_chromosome(self):
        pass

    def test_not_matching_reference(self):
        pass

    def test_wrong_transcript(self):
        pass

    def test_wrong_gene(self):
        pass

    def test_coordinate_overflow(self):
        pass

    def test_not_existing_patient(self):
        pass

    def test_batch_variants(self):
        pass

    def test_create_person(self):
        """
        Creates a person related to a patient
        :return:
        """
        # TODO: clarify why cannot send father and mother!
        person = {
            "patient_id": self.patient1_id,
            "relation": "brother",
            "relation_status": "unaffected"
        }
        person2 = {
            "patient_id": self.patient1_id,
            "relation": "sister",
            "relation_status": "unaffected"
        }
        persons = self.create_persons([person, person2])
        person_ids = [x["person_id"] for x in persons]
        for x in persons:
            self.assertTrue(x["person_id"] is not None, "Empty person id?")
            x.pop("person_id")
            self.assertTrue(x == person or x == person2, "Not equal to any person")
        self.delete_persons(person_ids)

    def test_create_phenotypes(self):
        """
        Creates a phenotype related to a patient
        :return:
        """
        person = {
            "patient_id": self.patient1_id,
            "relation": "brother",
            "relation_status": "unaffected"
        }
        persons = self.create_persons([person])
        person_ids = [x["person_id"] for x in persons]
        phenotype = {
            "person_id": person_ids[0],
            "phenotype_id": 4322,
            "observation": "present"
        }
        phenotypes = self.create_phenotypes([phenotype], person_ids[0])
        phenotype_ids = [x["person_phenotype_id"] for x in phenotypes]
        for x in phenotypes:
            self.assertTrue(x["person_phenotype_id"] is not None, "Empty person_phenotype id?")
            self.assertTrue(x["phenotype_id"] == phenotype["phenotype_id"], "Not equal to any phenotype")
            self.assertTrue(x["observation"] == phenotype["observation"], "Not equal to any phenotype")
        self.delete_phenotypes(phenotype_ids)
        self.delete_persons(person_ids)

    def create_phenotypes(self, phenotypes, person_id):
        response = requests.post("{url_base}persons/{person_id}/phenotypes".format(
            url_base=self.URL_BASE, person_id=person_id),
            headers=self.AUTH_HEADERS, data=json.dumps(phenotypes))
        self.assertTrue(response.status_code == 200, response.content)
        phenotype_ids = [x["person_phenotype_id"] for x in response.json()]
        output_phenotypes = []
        for phenotype_id in phenotype_ids:
            response = requests.get("{url_base}phenotypes/{phenotype_id}".format(
                url_base=self.URL_BASE, phenotype_id=phenotype_id), headers=self.AUTH_HEADERS)
            self.assertTrue(response.status_code == 200, response.content)
            output_phenotypes.append(response.json())
        return output_phenotypes

    def delete_phenotypes(self, phenotype_ids):
        for phenotype_id in phenotype_ids:
            response = requests.delete("{url_base}phenotypes/{phenotype_id}".format(
                url_base=self.URL_BASE, phenotype_id=phenotype_id),
                headers=self.AUTH_HEADERS)
            self.assertTrue(response.status_code == 200, response.content)

    def create_persons(self, persons):
        response = requests.post("{url_base}patients/{patient_id}/persons".format(
            url_base=self.URL_BASE, patient_id=self.patient1_id),
            headers=self.AUTH_HEADERS, data=json.dumps(persons))
        self.assertTrue(response.status_code == 200, response.content)
        person_ids = [x["person_id"] for x in response.json()]
        output_persons = []
        for person_id in person_ids:
            response = requests.get("{url_base}persons/{person_id}".format(
                url_base=self.URL_BASE, person_id=person_id), headers=self.AUTH_HEADERS)
            self.assertTrue(response.status_code == 200, response.content)
            output_persons.append(response.json())
        return output_persons

    def delete_persons(self, person_ids):
        for person_id in person_ids:
            response = requests.delete("{url_base}persons/{person_id}".format(
                url_base=self.URL_BASE, person_id=person_id),
                headers=self.AUTH_HEADERS)
            self.assertTrue(response.status_code == 200, response.content)

    def create_variants(self, variants, patient_id):
        """
        A list of variants that will be sent to Decipher
        :param variants:
        :param patient_id
        :return:
        """
        # send variants
        response = requests.post("{url_base}patients/{patient_id}/snvs".format(
            url_base=self.URL_BASE, patient_id=patient_id), headers=self.AUTH_HEADERS, data=json.dumps(variants))
        self.assertTrue(response.status_code == 200, response.content)
        variant_ids = [x["patient_snv_id"] for x in response.json()]

        # read variants
        output_variants = []
        for variant_id in variant_ids:
            response = requests.get("{url_base}snvs/{variant_id}".format(
                url_base=self.URL_BASE, variant_id=variant_id), headers=self.AUTH_HEADERS)
            self.assertTrue(response.status_code == 200, response.content)
            assert(response.status_code == 200)
            variant_json = response.json()
            output_variants.append(variant_json)
            assert(variant_json["patient_id"] == patient_id)
            assert(variant_json["assembly"] == variants[0]["assembly"])
            if len(variants) == 1:
                variant = variants[0]
                assert(variant_json["chr"] == variant["chr"])
                # NOTE: coordinates may change due to normalisation
                assert(variant_json["genotype"] == variant["genotype"])
                assert(variant_json["user_transcript"] == variant["user_transcript"])
                assert(variant_json["user_gene"] == variant["user_gene"])
                assert(variant_json["shared"] == variant["shared"])
        return output_variants

    def delete_variants(self, variant_ids):
        """
        Deletes variants from the database
        :param variant_ids:
        :return:
        """
        # delete variants
        for variant_id in variant_ids:
            response = requests.delete("{url_base}snvs/{variant_id}".format(
                url_base=self.URL_BASE, variant_id=variant_id), headers=self.AUTH_HEADERS)
            self.assertTrue(response.status_code == 200, response.content)
            assert(response.status_code == 200)

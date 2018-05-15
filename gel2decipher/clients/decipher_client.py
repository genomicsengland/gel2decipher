import logging
from gel2decipher.clients.rest_client import RestClient
from requests.exceptions import InvalidSchema


class DecipherClient(RestClient):

    def __init__(self, url_base, system_key, user_key):
        """
        System and user keys are required
        :param url_base:
        :param system_key:
        :param user_key:
        """
        RestClient.__init__(self, url_base)
        self.system_key = system_key
        self.user_key = user_key
        if not self.system_key or not self.user_key:
            raise ValueError("Authentication is required. Provide system and user key.")
        self.set_authenticated_header()
        # reads user and project
        response = self.get("info")
        self.project_id = response["user"]["project"]["project_id"]
        self.user_id = response["user"]["user_id"]
        logging.info(
            "Decipher client initialised with user id '{}' and project '{}'".format(self.user_id, self.project_id))

    def get_token(self):
        raise NotImplemented

    def set_authenticated_header(self, renew_token=False):
        self.headers["X-Auth-Token-System"] = self.system_key
        self.headers["X-Auth-Token-User"] = self.user_key

    def create_patients(self, patients):
        """

        :param patients: the list of input patients
        :type patients: list
        :return: the list of patient identifiers
        """
        if not all(patient.is_valid for patient in patients):
            validation_errors = [dict(patient.validation_errors) for patient in patients]
            raise InvalidSchema("Patients are invalid: {}".format(validation_errors), request=patients)
        response = self.post("projects/{project_id}/patients".format(project_id=self.project_id),
                             payload=[dict(patient) for patient in patients])
        return response

    def get_persons_by_patient(self, patient_id):
        response = self.get("patients/{patient_id}/persons".format(patient_id=patient_id))
        return response['persons']

    def delete_patient(self, patient_id):
        """

        :param patient_id:
        :return:
        """
        return self.delete("patients/{patient_id}".format(patient_id=patient_id))

    def create_persons(self, persons, patient_id):
        """

        :param persons:
        :type persons: list
        :param patient_id:
        :return:
        """
        if not all(person.is_valid for person in persons):
            validation_errors = [dict(patient.validation_errors) for patient in persons]
            raise InvalidSchema("Persons are invalid: {}".format(validation_errors), request=persons)
        response = self.post("patients/{patient_id}/persons".format(patient_id=patient_id),
                             payload=[dict(person) for person in persons])
        person_ids = [x["person_id"] for x in response]
        return person_ids

    def update_person(self, affection_status, person_id):
        """
        :type affection_status: str
        :type person_id: str
        :rtype:
        """
        response = self.patch("persons/{person_id}".format(person_id=person_id),
                              payload={"relation_status": affection_status})
        person_id = response["person_id"]
        return person_id

    def delete_person(self, person_id):
        """

        :param person_id:
        :return:
        """
        return self.delete("persons/{person_id}".format(person_id=person_id))

    def create_snvs(self, snvs, patient_id):
        """

        :param snvs:
        :param patient_id:
        :return:
        """
        if not all(snv.is_valid for snv in snvs):
            validation_errors = [dict(patient.validation_errors) for patient in snvs]
            raise InvalidSchema("Variants are invalid: {}".format(validation_errors), request=snvs)
        response = self.post("patients/{patient_id}/snvs".format(patient_id=patient_id),
                             payload=[dict(variant) for variant in snvs])
        variant_ids = [x["patient_snv_id"] for x in response]
        return variant_ids

    def get_snvs(self, patient_id):
        """

        :param patient_id:
        :return:
        """
        response = self.get("patients/{patient_id}/snvs".format(patient_id=patient_id))
        return response

    def delete_snv(self, snv_id):
        """

        :param snv_id:
        :return:
        """
        return self.delete("snvs/{snv_id}".format(snv_id=snv_id))

    def create_phenotypes(self, phenotypes, person_id):
        """

        :param phenotypes:
        :param person_id:
        :return:
        """
        if not all(phenotype.is_valid for phenotype in phenotypes):
            validation_errors = [dict(patient.validation_errors) for patient in phenotypes]
            raise InvalidSchema("Phenotypes are invalid: {}".format(validation_errors), request=phenotypes)
        response = self.post("persons/{person_id}/phenotypes".format(person_id=person_id),
                             payload=[dict(phenotype) for phenotype in phenotypes])
        phenotype_ids = [x["person_phenotype_id"] for x in response]
        return phenotype_ids

    def delete_phenotype(self, phenotype_id):
        """

        :param phenotype_id:
        :return:
        """
        return self.delete("phenotypes/{phenotype_id}".format(phenotype_id=phenotype_id))

# GEL 2 Decipher

GEL 2 Decipher allows you to send cases from Genomics England to Decipher by using several components:
* The CIPAPI REST API to fetch the pedigree
* The CVA REST API to fetch the variants
* The Decipher REST API to send all data over to Decipher system

## Creating persons from a pedigree

We are only creating the proband for any given family. 

There is some data we hold in pedigrees that is lost or transformed when mapping to Decipher's model:
* **Consent**. At Genomics England we have a consent based on 4 flags (programme consent, primary findings, secondary findings and carrier status). although only 2 are applicable to the proband. We are sending over only the consent for secondary findings.
* **Kariotypic sex**. At Genomics England we contemplate XXXY, XXX and XXYY, but this is a corner case as we have only have observed one case having XXXY.
* **Age**. We are inferring the age from the year of birth.
* **Inheritance**. Inheritance pattern cannot be sent over due to fundamental conflicts in the terms.

Technical debts:
* Create other members in the family
* Create phenotypes for other members of the family

## Sending phenotypes

No phenotypes were rejected. It is to be assessed if we need to run a migration of HPO terms.

In any case all HPO terms are added in the `note` field as free text.

## Sending variants

Technical debts:
* Indels are expected in VCF format?

## Selection of gene and transcript

1. The initial list of possible transcripts is read from Cellbase annotations in GRCh37
2. If we have a gene symbol provided by tiering we filter to only those transcripts belonging to that gene
3. Depending on the tier of the variant we select the transcripts for which our variant has a particular SO term:
    - Tier 1: "SO:0001893", "SO:0001574", "SO:0001575", "SO:0001587", "SO:0001589", "SO:0001578", "SO:0001582"
    - Tier 2: "SO:0001889", "SO:0001821", "SO:0001822", "SO:0001583", "SO:0001630", "SO:0001626"
4. We select transcripts having biotypes in: "IG_C_gene", "IG_D_gene"," IG_J_gene", "IG_V_gene", "IG_V_gene", 
    "protein_coding", "nonsense_mediated_decay", "non_stop_decay", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"
5. We select transcripts having a transcript flag: "basic"

The above procedure is equivalent to the procedure performed by tiering which means that we should not filter all transcripts in any case.

If after all filtering there are more than one transcript left we sort the Ensembl ids alphabetically and choose the first to ensure that our selection is deterministic.

We have a gene symbol provided by tiering process for each variant + all gene symbols of overlapping genes in GRCh37 provided by Cellbase

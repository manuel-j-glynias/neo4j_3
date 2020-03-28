import sys
import os
import json
import csv
import requests
from neo4j import GraphDatabase
import typing
import datetime

JAX_PATH: str = '/Users/mglynias/jax_march_2020/api-export/'

def send_mutation(mutation_payload:str, server:str) -> str:
    url = "http://" + server + ":7474/graphql/"
    headers = {
      'Authorization': 'Basic bmVvNGo6b21uaQ==',
      'Content-Type': 'application/json',
    }
    response = requests.request("POST", url, headers=headers, data = mutation_payload)
    if not response.ok:
        response.raise_for_status()
        sys.exit()

    responseBody: str = response.json()
    # print(responseBody)
    if 'errors' in responseBody:
        print(responseBody)
    return responseBody

def get_list_of_files(path: str) -> list:
    json_files = []
    for entry in os.scandir(path):
        if entry.is_file():
            json_files.append(entry.path)
    return json_files

def read_oncgenes_tumor_suppressors(file_name: str) -> dict:
    skip = True;
    gene_categories = {}
    with open(file_name) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            if skip:
                skip=False
            else:
                gene_categories[row[0]] = row[1]
    return gene_categories


def replace_characters(a_string: str):
    if a_string is not None:
        a_string = a_string.replace('α', 'alpha')
        a_string = a_string.replace('β', 'beta')
        a_string = a_string.replace('γ', 'gamma')
        a_string = a_string.replace('δ', 'delta')
        a_string = a_string.replace('ε', 'epsilon')
        a_string = a_string.replace('ζ', 'zeta')
        a_string = a_string.replace('η', 'eta')
        a_string = a_string.replace('θ', 'theta')
        a_string = a_string.replace('ι', 'iota')
        a_string = a_string.replace('κ', 'kappa')
        a_string = a_string.replace('σ', 'sigma')

        a_string = a_string.replace('ω', 'omega')
        a_string = a_string.replace(u"\u0394", 'delta')

        a_string = a_string.replace(u"\u025b", 'e')
        a_string = a_string.replace(u"\u223c", '~')

        a_string = a_string.replace("’", "")
        a_string = a_string.replace('"', '')
        a_string = a_string.replace(u"\u201c", '')
        a_string = a_string.replace(u"\u201d", '')
        a_string = a_string.replace(u"\u2018", '')
        a_string = a_string.replace(u"\u2019", '')
        a_string = a_string.replace(u"\u2020", '*')

        a_string = a_string.replace(u"\u207a", '+')
        a_string = a_string.replace(u"\u2011", '-')
        a_string = a_string.replace(u"\u2013", '-')
        a_string = a_string.replace(u"\u2014", '-')
        a_string = a_string.replace(u"\u2044", '/')
        a_string = a_string.replace(u"\u2122", 'TM')
        a_string = a_string.replace(u"\u2005", ' ')
        a_string = a_string.replace(u"\u2009", ' ')
        a_string = a_string.replace(u"\u2082", '2')
        a_string = a_string.replace('\t',' ')
        a_string = a_string.replace('\r','')
        a_string = a_string.replace('\n','')
        a_string = a_string.replace('⁸⁸','')

    return a_string

def read_one_gene_json(path,gene_categories) -> dict:
    gene_dict: dict = None
    with open(path, 'r') as afile:
        gene_data: dict = json.loads(afile.read())
        ckb_id: str = str(gene_data['id'])
        gene: str = str(gene_data['geneSymbol'])
        entrezId: str = str(gene_data['entrezId'])
        synonyms: list = gene_data['synonyms']
        chromosome: str = gene_data['chromosome']
        canonicalTranscript: str = gene_data['canonicalTranscript']
        if 'geneDescriptions' in gene_data:
            descriptions_array: list = gene_data['geneDescriptions']
            if len(descriptions_array)>0:
                description: str = replace_characters(descriptions_array[0]['description'])
                references: list = descriptions_array[0]['references']
                category = '-'
                if gene in gene_categories:
                    category = gene_categories[gene]
                gene_dict = {'ckb_id': ckb_id,
                        'gene':gene,
                        'description': description,
                        'references': references,
                        'entrezId' :entrezId,
                        'synonyms': synonyms,
                        'canonicalTranscript': canonicalTranscript,
                        'chromosome': chromosome,
                        'category':category}
    return gene_dict

# createJaxReference(
# PMID: String
# id: ID!
# jaxId: String!
# shortReference: String!
# title: String!): String

def write_one_reference(ref: dict,reference_dict:dict) -> (str,str):
    s = ''
    id = 'ref_' + str(ref['id'])
    title = ref['title']
    title = replace_characters(title)

    if not id in reference_dict:
        s = f'''{id} : createJaxReference(PMID: \\"{ref['pubMedId']}\\", id: \\"{id}\\", jaxId: \\"{str(ref['id'])}\\", shortReference: \\"{ref['url']}\\",title: \\"{title}\\"),'''
        reference_dict[id] = id
    return s, id





# createJaxGene(
# canonicalTranscript: [String]!
# chromosome: String!
# entrezId: String!
# id: ID!
# jaxId: String!
# name: String!
# statement: String!
# synonyms: [String]!): String

def write_one_gene(gene: dict,gene_dict:dict,reference_dict:dict)->str:
    id: str = 'g_' + gene['ckb_id']
    gene_dict[gene['gene']] = id
    synonyms: str = '['
    for syn in gene['synonyms']:
        synonyms += f'\\"{syn}\\",'
    synonyms += ']'

    s = f'''{id} : createJaxGene(canonicalTranscript: \\"{gene['canonicalTranscript']}\\", chromosome: \\"{gene['chromosome']}\\",entrezId: \\"{gene['entrezId']}\\", id: \\"{id}\\", jaxId: \\"{gene['ckb_id']}\\", name: \\"{gene['gene']}\\",oncogenic_category: \\"{gene['category']}\\", statement: \\"{gene['description']}\\",synonyms: {synonyms}),'''
    if len(gene['references'])>0:
        ref_array: str = '['
        for ref in gene['references']:
            m, ref_id = write_one_reference(ref,reference_dict)
            s += m
            ref_array += f'''\\"{ref_id}\\",'''
        ref_array += ']'
        #     addJaxGeneReferences(id: ID!references: [ID!]!): String

        s += f'''addJaxGeneReferences(id:\\"{id}\\", references:{ref_array}),'''
    return s


def write_genes(server:str,gene_dict:dict, reference_dict: dict) -> None:
    gene_categories = read_oncgenes_tumor_suppressors('data/OCP_oncogenes_tumorsuppressors.csv')
    json_files = get_list_of_files(JAX_PATH + 'genes')


    for json_file in json_files:
        gene: dict = read_one_gene_json(json_file, gene_categories)
        if gene is not None:
            mutation_payload: str = '{"query":"mutation {'
            mutation_payload += write_one_gene(gene,gene_dict,reference_dict)
            mutation_payload += '}"}'
            print(gene['gene'])
            send_mutation(mutation_payload,server)


def read_one_variant_json(path:str)->dict:
    with open(path, 'r') as afile:
        variant_data = json.loads(afile.read())
        ckb_id = str(variant_data['id'])
        full_name = variant_data['fullName']
        variant_type = variant_data['impact']
        protein_effect = variant_data['proteinEffect']
        description = '-'
        if 'geneVariantDescriptions' in variant_data:
            descriptions_array = variant_data['geneVariantDescriptions']
            if len(descriptions_array) > 0:
                description: str = replace_characters(descriptions_array[0]['description'])
                references: list = descriptions_array[0]['references']

        gene_id = variant_data['gene']['id']
        gene = variant_data['gene']['geneSymbol']
        gDot = '-'
        cDot = '-'
        pDot = '-'
        transcript = '-'
        if 'referenceTranscriptCoordinates' in variant_data and variant_data['referenceTranscriptCoordinates'] != None:
            gDot = variant_data['referenceTranscriptCoordinates']['gDna']
            cDot = variant_data['referenceTranscriptCoordinates']['cDna']
            pDot = variant_data['referenceTranscriptCoordinates']['protein']
            transcript = variant_data['referenceTranscriptCoordinates']['transcript']

        variant = {
            'ckb_id':ckb_id,
            'name':full_name,
            'variant_type':variant_type,
            'protein_effect':protein_effect,
            'description':description,
            'references': references,
            'gene_id':gene_id,
            'gene':gene,
            'gDot':gDot,
            'cDot':cDot,
            'pDot':pDot,
            'transcript':transcript,
        }
        return variant

def get_variant_id_from_jax_id(jax_id: str) -> str:
    return 'v_' + jax_id


# createJaxVariant(
# cDot: String!
# gDot: String!
# id: ID!
# jaxId: String!
# name: String!
# pDot: String!
# proteinEffect: String!
# statement: String!
# transcript: String!
# variant_type: String!): String

def write_one_variant(variant: dict,gene_dict:dict,reference_dict:dict) -> str:
    id: str = get_variant_id_from_jax_id(variant['ckb_id'])
    gene_id = gene_dict[variant['gene']]

    s = f'''{id} : createJaxVariant(cDot: \\"{variant['cDot']}\\",gDot: \\"{variant['gDot']}\\", id: \\"{id}\\", jaxId: \\"{variant['ckb_id']}\\", name: \\"{variant['name']}\\",pDot: \\"{variant['pDot']}\\",proteinEffect: \\"{variant['protein_effect']}\\", statement: \\"{variant['description']}\\", transcript: \\"{variant['transcript']}\\",variant_type: \\"{variant['variant_type']}\\"),'''
    # addJaxVariantGene(gene: [ID!]!id: ID!): String
    s += f'''addJaxVariantGene(gene: [\\"{gene_id}\\"], id:\\"{id}\\"),'''
    if len(variant['references']) > 0:
        ref_array: str = '['
        for ref in variant['references']:
            m, ref_id = write_one_reference(ref, reference_dict)
            s += m
            ref_array += f'''\\"{ref_id}\\",'''
        ref_array += ']'
        #     addJaxGeneReferences(id: ID!references: [ID!]!): String

        s += f'''addJaxVariantReferences(id:\\"{id}\\", references:{ref_array}),'''
    return s


def write_variants(server:str,gene_dict:dict,reference_dict:dict) -> None:
    json_files = get_list_of_files(JAX_PATH + 'variants')
    for json_file in json_files:
        variant: dict = read_one_variant_json(json_file)
        if variant is not None:
            mutation_payload: str = '{"query":"mutation {'
            mutation_payload += write_one_variant(variant,gene_dict,reference_dict)
            mutation_payload += '}"}'
            print(variant['name'],variant['ckb_id'])
            if variant['name']=='BAP1 wild-type':
                print(variant)
            send_mutation(mutation_payload,server)

def read_disease_ontology(filepath:str):
    ontology_dict = {}
    ontology_names_dict = {}
    children_dict = {}
    is_a = []
    doid = ''
    with open(filepath) as fp:
        for line in fp:
            if line.startswith('[Term]'):
                if len(is_a) > 0:
                    ontology_dict[doid] = is_a
                is_a = []
            elif line.startswith('id'):
                # id: DOID:0001816
                doid = line[4:].rstrip()
            elif line.startswith('name'):
                name = line[6:].rstrip()
                ontology_names_dict[doid] = name
            elif line.startswith('is_a'):
                # is_a: DOID:175 ! vascular cancer
                parent = line[6:line.find('!') - 1]
                is_a.append(parent)
                if parent in children_dict:
                    children = children_dict[parent]
                    children.append(doid)
                else:
                    children_dict[parent] = [doid]

    return ontology_dict, ontology_names_dict, children_dict


def show_nodes(id:str,ontology_dict, ontology_names_dict, children_dict,depth):
    name = id
    if id in ontology_names_dict:
        name = ontology_names_dict[id]
    print('-' * depth,name)
    if id in children_dict:
        for node in children_dict[id]:
            show_nodes(node,ontology_dict, ontology_names_dict, children_dict,depth+1)

# createJaxDisease(id: ID!name: String!): String

# addJaxDiseaseChildren(children: [ID!]!id: ID!): String

def get_ID_from_doid(doid:str)->str:
    return 'disease_' + doid[5:]

def get_disease_id_from_jax_id(doid:str):
    if doid == 'JAX:10000003':  # substitute id for 'cancer' for 'Advanced Solid Tumor" to get children
        doid = 'DOID:162'
    return get_ID_from_doid(doid)

def write_nodes(doid:str,ontology_dict, ontology_names_dict, children_dict,server,id_dict)->None:
    name = doid
    if doid in ontology_names_dict:
        name = ontology_names_dict[doid]
    id = get_ID_from_doid(doid)
    children = '['
    if doid in children_dict:
        for node in children_dict[doid]:
            write_nodes(node,ontology_dict, ontology_names_dict, children_dict,server,id_dict)
            child_id = get_ID_from_doid(node)
            children += f'''\\"{child_id}\\"'''
    children += ']'

    mutation_payload: str = ''
    if id not in id_dict:
        mutation_payload += f'''{id}: createJaxDisease(id:\\"{id}\\", name:\\"{name}\\"),'''
        id_dict[id] = id
    if doid in children_dict:
        mutation_payload += f'''addJaxDiseaseChildren(children:{children}, id:\\"{id}\\"),'''
    if len(mutation_payload)>0:
        mutation_payload = '{"query":"mutation {' + mutation_payload + '}"}'
        send_mutation(mutation_payload, server)


def write_diseases(server:str) -> None:
    ontology_dict, ontology_names_dict, children_dict = read_disease_ontology('data/doid.obo')
    id_dict: dict = {}
    write_nodes('DOID:162',ontology_dict, ontology_names_dict, children_dict,server,id_dict)

def get_drug_class_id_from_jax_id(jax_id:str)->str:
    return 'drugclass_' + jax_id


# createJaxDrugClass(
# id: ID!
# jaxId: String!
# name: String!): String

def read_one_drug_class_json(path) -> dict:
    drug_class: dict = None
    with open(path, 'r') as afile:
        drug_class_data = json.loads(afile.read())
        id: str = get_drug_class_id_from_jax_id(str(drug_class_data['id']))
        name: str = str(drug_class_data['drugClass'])
        jaxId: str = str(drug_class_data['id'])
        drug_class = {'id': id, 'name':name, 'jaxId':jaxId}
    return drug_class

# createJaxDrugClass(
# id: ID!
# jaxId: String!
# name: String!): String

def write_drug_classes(server:str) -> None:
    json_files = get_list_of_files(JAX_PATH + 'drugClasses')
    for json_file in json_files:
        drug_class: dict = read_one_drug_class_json(json_file)
        if drug_class is not None:
            mutation_payload: str = '{"query":"mutation {'
            mutation_payload += f'''{drug_class['id']}: createJaxDrugClass(id:\\"{drug_class['id']}\\",jaxId:\\"{drug_class['jaxId']}\\", name:\\"{drug_class['name']}\\"),'''
            mutation_payload += '}"}'
            send_mutation(mutation_payload,server)


def get_drug_id_from_jax_id(jax_id:str)->str:
    return 'drug_' + jax_id


def read_one_drug_json(path) -> dict:
    drug_dict: dict = None
    with open(path, 'r') as afile:
        drug_data: dict = json.loads(afile.read())
        jax_id: str = str(drug_data['id'])
        id: str = get_drug_id_from_jax_id(jax_id)
        name: str = str(drug_data['drugName'])
        synonyms: list = drug_data['synonyms']
        tradeName: str = drug_data['tradeName']
        drugClasses: list = drug_data['drugClasses']
        if 'drugDescriptions' in drug_data:
            descriptions_array: list = drug_data['drugDescriptions']
            if len(descriptions_array)>0:
                description: str = replace_characters(descriptions_array[0]['description'])
                references: list = descriptions_array[0]['references']
                drug_dict = {'jax_id': jax_id,
                        'description': description,
                        'references': references,
                        'id' :id,
                        'name':name,
                        'synonyms': synonyms,
                        'tradeName': tradeName,
                        'drugClasses': drugClasses,
                        }
    return drug_dict


# createJaxDrug(
# id: ID!
# jaxId: String!
# name: String!
# statement: String!
# synonyms: [String]!
# tradeName: String!): String

def write_one_drug(drug_dict: dict,reference_dict: dict) -> str:
    synonyms: str = '['
    for syn in drug_dict['synonyms']:
        synonyms += f'\\"{syn}\\",'
    synonyms += ']'
    s = f'''{drug_dict['id']} : createJaxDrug(id: \\"{drug_dict['id']}\\", jaxId: \\"{drug_dict['jax_id']}\\", name: \\"{drug_dict['name']}\\",statement: \\"{drug_dict['description']}\\", synonyms: {synonyms},tradeName: \\"{drug_dict['tradeName']}\\"),'''

    # addJaxDrugDrugClasses(drugClasses: [ID!]!id: ID!): String
    # Adds DrugClasses to JaxDrug entity
    drugClasses: str = '['
    for dc in drug_dict['drugClasses']:
        id:str = get_drug_class_id_from_jax_id(str(dc['id']))
        drugClasses += f'\\"{id}\\",'
    drugClasses += ']'

    s += f'''addJaxDrugDrugClasses(drugClasses: {drugClasses}, id:\\"{drug_dict['id']}\\"),'''

    if len(drug_dict['references']) > 0:
        ref_array: str = '['
        for ref in drug_dict['references']:
            m, ref_id = write_one_reference(ref, reference_dict)
            s += m
            ref_array += f'''\\"{ref_id}\\",'''
        ref_array += ']'
        #     addJaxDrugReferences(id: ID!references: [ID!]!): String

        s += f'''addJaxDrugReferences(id:\\"{drug_dict['id']}\\", references:{ref_array}),'''
    return s


def write_drugs(server:str,reference_dict:dict) -> None:
    json_files = get_list_of_files(JAX_PATH + 'drugs')
    for json_file in json_files:
        drug_dict: dict = read_one_drug_json(json_file)
        if drug_dict is not None:
            mutation_payload: str = '{"query":"mutation {'
            mutation_payload += write_one_drug(drug_dict,reference_dict)
            mutation_payload += '}"}'
            print(drug_dict['name'])
            send_mutation(mutation_payload,server)


def get_therapy_id_from_jax_id(jax_id:str):
    return 'therapy_' + jax_id


def read_one_therapy_json(path:str):
    therapy_dict: dict = None
    with open(path, 'r') as afile:
        therapy_data = json.loads(afile.read())
        id: str = get_therapy_id_from_jax_id(str(therapy_data['id']))
        name: str = str(therapy_data['therapyName'])
        jaxId: str = str(therapy_data['id'])
        drugs: list = therapy_data['drugs']
        therapy_dict = {'id': id, 'name': name, 'jaxId': jaxId, 'drugs':drugs}
    return therapy_dict


# createJaxTherapy(
# id: ID!
# jaxId: String!
# name: String!): String
# Creates a JaxTherapy entity

def write_therapies(server:str) -> None:
    json_files = get_list_of_files(JAX_PATH + 'therapies')
    for json_file in json_files:
        therapy_dict: dict = read_one_therapy_json(json_file)
        if therapy_dict is not None:
            mutation_payload: str = '{"query":"mutation {'
            mutation_payload += f'''{therapy_dict['id']}: createJaxTherapy(id:\\"{therapy_dict['id']}\\",jaxId:\\"{therapy_dict['jaxId']}\\", name:\\"{therapy_dict['name']}\\"),'''

            # addJaxTherapyDrugs(drugs: [ID!]!id: ID!): String
            # Adds Drugs to JaxTherapy entity
            drugs: str = '['
            for dc in therapy_dict['drugs']:
                id: str = get_drug_id_from_jax_id(str(dc['id']))
                drugs += f'\\"{id}\\",'
            drugs += ']'
            mutation_payload += f'''addJaxTherapyDrugs(drugs: {drugs}, id:\\"{therapy_dict['id']}\\"),'''

            mutation_payload += '}"}'
            print(therapy_dict['name'])
            send_mutation(mutation_payload, server)


def get_molecular_profile_id_from_jax_id(jax_id:str):
    return 'mp_' + jax_id

def read_one_molecular_profile_json(path:str):
    mp_dict: dict = None
    with open(path, 'r') as afile:
        mp_data = json.loads(afile.read())
        id: str = get_molecular_profile_id_from_jax_id(str(mp_data['id']))
        name: str = replace_characters(str(mp_data['profileName']))
        jaxId: str = str(mp_data['id'])
        variants: list = mp_data['geneVariants']
        mp_dict = {'id': id, 'name': name, 'jaxId': jaxId, 'variants':variants}
    return mp_dict

# createJaxMolecularProfile(
# id: ID!
# jaxId: String!
# name: String!): String
# Creates a JaxMolecularProfile entity

def write_molecular_profiles(server:str) -> None:
    json_files = get_list_of_files(JAX_PATH + 'molecularProfiles')
    for json_file in json_files:
        mp_dict: dict = read_one_molecular_profile_json(json_file)
        if mp_dict is not None:
            mutation_payload: str = '{"query":"mutation {'
            mutation_payload += f'''{mp_dict['id']}: createJaxMolecularProfile(id:\\"{mp_dict['id']}\\",jaxId:\\"{mp_dict['jaxId']}\\", name:\\"{mp_dict['name']}\\"),'''
             # addJaxMolecularProfileVariants(id: ID!variants: [ID!]!): String
            # Adds Variants to JaxMolecularProfile entity
            variants: str = '['
            for var in mp_dict['variants']:
                id: str = get_variant_id_from_jax_id(str(var['id']))
                variants += f'\\"{id}\\",'
            variants += ']'
            mutation_payload += f'''addJaxMolecularProfileVariants(id:\\"{mp_dict['id']}\\", variants: {variants} ),'''

            mutation_payload += '}"}'
            print(mp_dict['name'])
            if mp_dict['name'] == 'ARID1A E1776* ARID1A 	S705fs ESR1 D538G PIK3CA N345K':
                print(mp_dict)
            send_mutation(mutation_payload, server)


status = {'FDA approved': 1,
          'FDA approved - On Companion Diagnostic': 1,
          'FDA approved - Has Companion Diagnostic': 1,
          'FDA contraindicated': 2,
          'Guideline': 3,
          'Phase III': 4,
          'Phase II': 5,
          'Phase Ib/II': 5,
          'Phase I': 6,
          'Phase 0': 6,
          'Case Reports/Case Series': 7,
          'Clinical Study - Meta-analysis': 7,
          'Clinical Study - Cohort': 7,
          'Clinical Study': 7,
          'Preclinical': 8,
          'Preclinical - Patient cell culture': 8,
          'Preclinical - Pdx': 8,
          'Preclinical - Pdx & cell culture': 8,
          'Preclinical - Cell culture': 8,
          'Preclinical - Cell line xenograft': 8
          }


def get_approval_level_from_status(status_string:str)->int:
    level: int = 9
    if status_string in status:
        level = status[status_string]
    return level


def read_one_indication_json(path:str):
    indication_dict_list: list = []
    with open(path, 'r') as afile:
        indication_data = json.loads(afile.read())
        disease_id: str = get_disease_id_from_jax_id(indication_data['termId'])
        disease_name: str = indication_data['name']
        for result in indication_data['evidence']:
            indication_dict: dict = {
                'disease_id': disease_id,
                'molecular_profile_id': get_molecular_profile_id_from_jax_id(str(result['molecularProfile']['id'])),
                'therapy_id': get_therapy_id_from_jax_id(str(result['therapy']['id'])),
                'ampCapAscoEvidenceLevel': result['ampCapAscoEvidenceLevel'],
                'ampCapAscoInferredTier': result['ampCapAscoInferredTier'],
                'approvalLevel':  get_approval_level_from_status(result['approvalStatus']),
                'approvalStatus': result['approvalStatus'],
                'evidenceType': result['evidenceType'],
                'id': 'result_' + str(result['id']),
                'name': replace_characters(disease_name + '|' + result['molecularProfile']['profileName'] + '|' + result['therapy']['therapyName']),
                'responseType': result['responseType'],
                'statement': replace_characters(result['efficacyEvidence']),
                'references': result['references']
            }
            indication_dict_list.append(indication_dict)
    return indication_dict_list
# createJaxResult(
# ampCapAscoEvidenceLevel: String!
# ampCapAscoInferredTier: String!
# approvalLevel: Int!
# approvalStatus: String!
# evidenceType: String!
# id: ID!
# name: String!
# responseType: String!
# statement: String!): String
# Creates a JaxResult entity

# addJaxResultDisease(disease: [ID!]!id: ID!): String
# Adds Disease to JaxResult entity
#
# addJaxResultMolecular_profile(id: ID!molecular_profile: [ID!]!): String
# Adds Molecular_profile to JaxResult entity
#
# addJaxResultReferences(id: ID!references: [ID!]!): String
# Adds References to JaxResult entity
#
# addJaxResultTherapy(id: ID!therapy: [ID!]!): String
# Adds Therapy to JaxResult entity

def write_results(server:str, reference_dict:dict) -> None:
    json_files = get_list_of_files(JAX_PATH + 'indications')
    count: int = 0
    for json_file in json_files:
        indication_dict_list: list = read_one_indication_json(json_file)
        for result in indication_dict_list:
            mutation_payload: str = '{"query":"mutation {'
            mutation_payload += f'''{result['id']}: createJaxResult(ampCapAscoEvidenceLevel:\\"{result['ampCapAscoEvidenceLevel']}\\",ampCapAscoInferredTier:\\"{result['ampCapAscoInferredTier']}\\",approvalLevel:{result['approvalLevel']},approvalStatus:\\"{result['approvalStatus']}\\",evidenceType:\\"{result['evidenceType']}\\",id:\\"{result['id']}\\",name:\\"{result['name']}\\",responseType:\\"{result['responseType']}\\",statement:\\"{result['statement']}\\",),'''
            mutation_payload += f'''addJaxResultDisease(disease: [\\"{result['disease_id']}\\"],id:\\"{result['id']}\\" ),'''
            mutation_payload += f'''addJaxResultMolecular_profile(id:\\"{result['id']}\\", molecular_profile: [\\"{result['molecular_profile_id']}\\"], ),'''
            mutation_payload += f'''addJaxResultTherapy(id:\\"{result['id']}\\", therapy: [\\"{result['therapy_id']}\\"], ),'''
            if len(result['references']) > 0:
                ref_array: str = '['
                for ref in result['references']:
                    m, ref_id = write_one_reference(ref, reference_dict)
                    mutation_payload += m
                    ref_array += f'''\\"{ref_id}\\",'''
                ref_array += ']'
                mutation_payload += f'''addJaxResultReferences(id:\\"{result['id']}\\", references:{ref_array}),'''

            mutation_payload += '}"}'
            print(count, result['name'])
            count = count + 1
            # if result['name'] == 'ARID1A E1776* ARID1A 	S705fs ESR1 D538G PIK3CA N345K':
            #     print(result)
            send_mutation(mutation_payload, server)

def erase_neo4j(server):
    uri = "bolt://" + server + ":7687"
    with open('schema.graphql', 'r') as file:
        idl_as_string = file.read()
    driver = GraphDatabase.driver(uri, auth=("neo4j", "omni"))
    with driver.session() as session:
        tx = session.begin_transaction()
        tx.run("match(a) detach delete(a)")
        result = tx.run("call graphql.idl('" + idl_as_string + "')")
        print(result.single()[0])
        tx.commit()
    driver.close()

def main():
    server:str = 'localhost'
    # server: str = '165.227.89.140'
    erase_neo4j(server)

    gene_dict: dict = {}
    reference_dict: dict = {}

    write_genes(server,gene_dict,reference_dict)
    write_variants(server,gene_dict,reference_dict)
    write_diseases(server)
    write_drug_classes(server)
    write_drugs(server, reference_dict)
    write_therapies(server)
    write_molecular_profiles(server)
    write_results(server, reference_dict)


if __name__ == "__main__":
    main()
#
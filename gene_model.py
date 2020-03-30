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

def read_one_gene_jax_json(path) -> dict:
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
                gene_dict = {'ckb_id': ckb_id,
                        'gene':gene,
                        'description': description,
                        'references': references,
                        'entrezId' :entrezId,
                        'synonyms': synonyms,
                        'canonicalTranscript': canonicalTranscript,
                        'chromosome': chromosome}
    return gene_dict

# createJaxReference(
# PMID: String
# id: ID!
# jaxId: String!
# shortReference: String!
# title: String!): String

def write_one_jax_reference(ref: dict,reference_dict:dict) -> (str,str):
    s = ''
    id = 'ref_' + str(ref['id'])
    title = ref['title']
    title = replace_characters(title)

    if not id in reference_dict:
        s = f'''{id} : createJaxReference(PMID: \\"{ref['pubMedId']}\\", id: \\"{id}\\", jaxId: \\"{str(ref['id'])}\\", shortReference: \\"{ref['url']}\\",title: \\"{title}\\"),'''
        reference_dict[id] = id
    return s, id


def get_gene_id_from_jax_id(jax_id:str)->str:
    return 'jax_gene_' +jax_id


# createJaxGene(
# canonicalTranscript: [String]!
# chromosome: String!
# entrezId: String!
# id: ID!
# jaxId: String!
# name: String!
# statement: String!
# synonyms: [String]!): String

def write_one_jax_gene(gene: dict,gene_dict:dict,reference_dict:dict)->str:
    id: str = get_gene_id_from_jax_id(str(gene['ckb_id']))
    gene_dict[gene['gene']] = id
    synonyms: str = '['
    for syn in gene['synonyms']:
        synonyms += f'\\"{syn}\\",'
    synonyms += ']'

    s = f'''{id} : createJaxGene(canonicalTranscript: \\"{gene['canonicalTranscript']}\\", chromosome: \\"{gene['chromosome']}\\",entrezId: \\"{gene['entrezId']}\\", id: \\"{id}\\", jaxId: \\"{gene['ckb_id']}\\", name: \\"{gene['gene']}\\", statement: \\"{gene['description']}\\",synonyms: {synonyms}),'''
    if len(gene['references'])>0:
        ref_array: str = '['
        for ref in gene['references']:
            m, ref_id = write_one_jax_reference(ref,reference_dict)
            s += m
            ref_array += f'''\\"{ref_id}\\",'''
        ref_array += ']'
        #     addJaxGeneReferences(id: ID!references: [ID!]!): String

        s += f'''addJaxGeneReferences(id:\\"{id}\\", references:{ref_array}),'''
    return s


def write_jax_genes(server:str,gene_dict:dict, reference_dict: dict) -> None:
    # gene_categories = read_oncgenes_tumor_suppressors('data/OCP_oncogenes_tumorsuppressors.csv')
    json_files = get_list_of_files(JAX_PATH + 'genes')

    for json_file in json_files:
        gene: dict = read_one_gene_jax_json(json_file)
        if gene is not None:
            mutation_payload: str = '{"query":"mutation {'
            mutation_payload += write_one_jax_gene(gene,gene_dict,reference_dict)
            mutation_payload += '}"}'
            print(gene['gene'])
            send_mutation(mutation_payload,server)


def read_omni_genes(file_name: str)->list:
    gene_categories = []
    with open(file_name) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene_data:dict = {'symbol':row['HGNC_Symbol'],
                              'category':row['GeneType'],
                              'category_source': row['GeneTypeSource']}

            gene_categories.append(gene_data)
    return gene_categories


def fetch_gene_id_by_gene_name(gene_name:str)->str:
    entrezgene_id = None
    requestURL = "http://mygene.info/v3/query?species=human&q=symbol:" + gene_name
    r = requests.get(requestURL, headers={"Accept": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = r.json()
    # print(responseBody)
    if 'hits' in responseBody and len(responseBody['hits'])>0 and 'entrezgene' in responseBody['hits'][0]:
        entrezgene_id = responseBody['hits'][0]['entrezgene']
    return entrezgene_id


def fetch_gene_info_by_gene_id(gene_id:str)->dict:
  requestURL = "http://mygene.info/v3/gene/" + gene_id
  r = requests.get(requestURL, headers={ "Accept" : "application/json"})
  if not r.ok:
    r.raise_for_status()
    sys.exit()

  responseBody = r.json()
  # pprint.pprint(responseBody)
  return responseBody


def erase_neo4j(server):
    uri = "bolt://" + server + ":7687"
    with open('schema2.graphql', 'r') as file:
        idl_as_string = file.read()
    driver = GraphDatabase.driver(uri, auth=("neo4j", "omni"))
    with driver.session() as session:
        tx = session.begin_transaction()
        tx.run("match(a) detach delete(a)")
        result = tx.run("call graphql.idl('" + idl_as_string + "')")
        print(result.single()[0])
        tx.commit()
    driver.close()

def create_mygene_reference(gene_id):
  reference = {}
  reference['type'] = 'InternetReference'
  reference['url'] = "http://mygene.info/v3/gene/" + gene_id
  reference['accessed_date'] = datetime.datetime.now()
  return reference

def create_hgnc_gene_name_dict()->dict:
    genes = {}
    genes['C10orf54'] = 'VSIR'
    genes['C11orf30'] = 'EMSY'
    genes['FAM123B'] = 'AMER1'
    genes['FAM175A'] = 'ABRAXAS1'
    genes['FAM46C'] = 'TENT5C'
    genes['GPR124'] = 'ADGRA2'
    genes['H3F3A'] = 'H3-3A'
    genes['H3F3B'] = 'H3-3B'
    genes['H3F3C'] = 'H3-5'
    genes['HIST1H1C'] = 'H1-2'
    genes['HIST1H3B'] = 'H3C2'
    genes['HIST1H2BB'] = 'H2BC3'
    genes['HIST1H2BD'] = 'H2BC5'
    genes['HIST1H3A'] = 'H3C1'
    genes['HIST1H3E'] = 'H3C6'
    genes['HIST1H3F'] = 'H3C7'
    genes['HIST1H3G'] = 'H3C8'
    genes['HIST1H3H'] = 'H3C10'
    genes['HIST1H3I'] = 'H3C11'
    genes['HIST1H3J'] = 'H3C12'
    genes['HIST2H3A'] = 'H3C15'
    genes['HIST2H3C'] = 'H3C14'
    genes['HIST2H3D'] = 'H3C13'
    genes['HIST3H3'] = 'H3-4'
    genes['LOC101926927'] = 'withdrawn'
    genes['MARCH9'] = 'MARCHF9'
    genes['MLL'] = 'KMT2A'
    genes['MRE11A'] = 'MRE11'
    genes['MYCL1'] = 'MYCL'
    genes['PAK7'] = 'PAK5'
    genes['PARK2'] = 'PRKN'
    genes['RFWD2'] = 'COP1'
    genes['TCEB1'] = 'ELOC'
    genes['WISP3'] = 'CCN6'
    return genes

def get_name_for_internet_reference(url:str ,accessed_date: str):
    pos1 = url.find('//') + 2
    pos2 = url.find('/',pos1)
    name = url[pos1:pos2]
    name += ' (accessed on:' + accessed_date + ')'
    return name


def get_acessed_date_as_string(d:datetime)-> str:
    date_time = d.strftime("%m/%d/%Y")
    return date_time


# createMyGeneInfo_Gene(
# chromosome: String!
# end: Int
# entrezId: String!
# id: ID!
# name: String!
# start: Int
# statement: String!
# strand: Strand!
# synonyms: [String]!): String
# Creates a MyGeneInfo_Gene entity


def get_gene_id_from_entrez_id(entrez_id:str)->str:
    return 'geneInfo_gene_' + entrez_id


def create_myGeneInfo_gene(omni_gene:dict, server:str)->None:
    id = get_gene_id_from_entrez_id(omni_gene['entrez_gene_id'])
    gene: str = omni_gene['symbol']
    chrom: str = omni_gene['chrom']
    strand: str = omni_gene['strand']
    start: int = omni_gene['start']
    end: int = omni_gene['end']
    entrez_id = omni_gene['entrez_gene_id']
    statement: str = gene + ' on chromosome ' + chrom + ' at ' + str(start) + '-' + str(end)
    if 'summary' in omni_gene:
        statement = omni_gene['summary']
    synonyms: str = '['
    for syn in omni_gene['synonyms']:
        synonyms += f'\\"{syn}\\",'
    synonyms += ']'
    mutation_payload: str = '{"query":"mutation {'
    s = f'{id}: createMyGeneInfo_Gene(chromosome: \\"{chrom}\\", end: {end}, entrezId: \\"{entrez_id}\\", id: \\"{id}\\", name: \\"{gene}\\" start: {start}, statement:\\"{statement}\\", strand:{strand},synonyms: {synonyms}),'
    ref = omni_gene['reference']
    if ref['type'] == 'InternetReference':
        ref_id: str = 'ref_' + id
        accessed: str = get_acessed_date_as_string(ref['accessed_date'])
        ir_name: str = get_name_for_internet_reference(ref["url"], accessed)
        # createInternetReference(
        # accessed_date: String!
        # id: ID!
        # shortReference: String!
        # web_address: String!): String
        # Creates a InternetReference entity
        s += f'{ref_id}: createInternetReference(accessed_date:\\"{accessed}\\", id:\\"{ref_id}\\", shortReference: \\"{ir_name}\\", web_address:\\"{ref["url"]}\\" ),'
        ref_id2 = 'gref_' + id
        s += f'{ref_id2}: addMyGeneInfo_GeneReferences(id:\\"{id}\\", references:[\\"{ref_id}\\"] ),'
    mutation_payload += s + '}"}'
    print(mutation_payload)
    send_mutation(mutation_payload, server)


def get_omnigene_id_from_entrez_id(entrez_id:str)->str:
    return 'omnigene_' + entrez_id

# createEditableStatement(
# deleted: Boolean!
# edit_date: String!
# editor: String!
# field: String!
# id: ID!
# statement: String!): String

def createEditableStatement(statement:str, field:str, editor:str, esList:list) -> (str,str):
    now = datetime.datetime.now()
    edit_date:str = now.strftime("%Y-%m-%d-%H-%M-%S-%f")
    numES: int = len(esList)
    id:str = 'es_' + str(numES)
    esList.append(id)
    s = f'''{id} : createEditableStatement(deleted: false, edit_date: \\"{edit_date}\\", editor: \\"{editor}\\",field: \\"{field}\\", id: \\"{id}\\",statement: \\"{statement}\\"),'''
    return s, id


# createOmniGene(
# id: ID!
# name: String!
# panelName: String!): String
# Creates a OmniGene entity

def create_omniGene(omni_gene:dict, gene_dict:dict,esList, editor, server:str)->None:
    id = get_omnigene_id_from_entrez_id(omni_gene['entrez_gene_id'])
    gene: str = omni_gene['symbol']
    panel_name = omni_gene['panel_name']
    mutation_payload: str = '{"query":"mutation {'
    s = f'{id}: createOmniGene(id: \\"{id}\\", name: \\"{gene}\\", panelName:\\"{panel_name}\\" ),'
    # create geneDescription EditableStatement
    field1: str = 'geneDescription_' + id
    statement1: str = '(Insert Gene Description)'
    (m, id1) = createEditableStatement(statement1,field1,editor,esList)
    s += m
#     addOmniGeneGeneDescription(geneDescription: [ID!]!id: ID!): String
#     Adds GeneDescription to OmniGene entity
    s += f'addOmniGeneGeneDescription(geneDescription:[\\"{id1}\\"], id:\\"{id}\\" ),'
    # create OncogenicCategory EditableStatement
    field2: str = 'OncogenicCategory_' + id
    statement2: str = '(Draft) ' + omni_gene['category']
    (m, id2) = createEditableStatement(statement2,field2,editor,esList)
    s += m
 # addOmniGeneOncogenicCategory(id: ID!oncogenicCategory: [ID!]!): String
# Adds OncogenicCategory to OmniGene entity
    s += f'addOmniGeneOncogenicCategory(id:\\"{id}\\", oncogenicCategory:[\\"{id2}\\"] ),'
    # addOmniGeneSynonymsString(id: ID!synonymsString: [ID!]!): String
    # Adds SynonymsString to OmniGene entity
    field3: str = 'SynonymsString_' + id
    statement3: str = ''
    if 'synonym' in omni_gene:
        statement3 = omni_gene['synonym']
    (m, id3) = createEditableStatement(statement3, field3, editor, esList)
    s += m
    s += f'addOmniGeneSynonymsString(id:\\"{id}\\", synonymsString:[\\"{id3}\\"] ),'

    # addOmniGeneJaxGene(id: ID!jaxGene: [ID!]!): String
# Adds JaxGene to OmniGene entity
#     jaxGene = get_gene_id_from_jax_id(omni_gene['entrez_gene_id'])
    if gene in gene_dict:
        jaxGene = gene_dict[gene]
        s += f'addOmniGeneJaxGene(id:\\"{id}\\", jaxGene:[\\"{jaxGene}\\"] ),'
    else:
        print("no jax gene for ",gene)
# addOmniGeneMyGeneInfoGene(id: ID!myGeneInfoGene: [ID!]!): String
# Adds MyGeneInfoGene to OmniGene entity
    myGeneInfoGene = get_gene_id_from_entrez_id(omni_gene['entrez_gene_id'])
    s += f'addOmniGeneMyGeneInfoGene(id:\\"{id}\\", myGeneInfoGene:[\\"{myGeneInfoGene}\\"] ),'

    mutation_payload: str = '{"query":"mutation {'
    mutation_payload += s + '}"}'
    # print(mutation_payload)
    send_mutation(mutation_payload, server)


def write_omni_genes(omini_genes:list, gene_dict:dict, server:str):
    hgnc_gene_name_dict = create_hgnc_gene_name_dict()
    esList: list = []
    for omni_gene in omini_genes:
        gene_name: str = omni_gene['symbol']
        omni_gene['panel_name'] = gene_name
        if gene_name in hgnc_gene_name_dict:
            omni_gene['panel_name'] = gene_name
            omni_gene['synonym'] = gene_name
            gene_name = hgnc_gene_name_dict[gene_name]
            omni_gene['symbol'] = gene_name
        entrez_gene_id = fetch_gene_id_by_gene_name(gene_name)
        omni_gene['entrez_gene_id'] = entrez_gene_id
        if entrez_gene_id is None:
            print("no entrz gene id for",gene_name)
        else:
            gene_info = fetch_gene_info_by_gene_id(entrez_gene_id)
            if 'genomic_pos_hg19' in gene_info:
                genomic_pos = gene_info['genomic_pos_hg19']
                if isinstance(genomic_pos, list):
                    genomic_pos = genomic_pos[0]
                omni_gene['chrom'] = genomic_pos['chr']
                omni_gene['strand'] = 'FORWARD'
                if genomic_pos['strand'] == -1:
                    omni_gene['strand'] = 'REVERSE'
                omni_gene['start'] = genomic_pos['start']
                omni_gene['end'] = genomic_pos['end']
                omni_gene['reference'] = create_mygene_reference(entrez_gene_id)
                if 'summary' in gene_info:
                    omni_gene['summary'] = gene_info['summary']
                if 'alias' in gene_info:
                    omni_gene['synonyms'] = gene_info['alias']
                omni_gene['synonyms'] = []
                print(gene_name)
                create_myGeneInfo_gene(omni_gene,server)
                create_omniGene(omni_gene,gene_dict,esList,'auto',server)
            else:
                print('no gene_info for',gene_name)


def main():
    server:str = 'localhost'
    # server: str = '165.227.89.140'
    erase_neo4j(server)

    gene_dict: dict = {}
    reference_dict: dict = {}

    write_jax_genes(server,gene_dict,reference_dict)
    omini_genes:list = read_omni_genes('data/tblOS_GLOBAL_GLOBAL_Ref_AllPanelsGenes.csv')
    write_omni_genes(omini_genes,gene_dict,server)



if __name__ == "__main__":
    main()
#
import csv
import requests
from neo4j import GraphDatabase
import datetime
import sys
import os
os.environ['NCBI_API_KEY'] = 'cde5c1a63fa16711994bfe74b858747cbb08'
from metapub import PubMedFetcher
import re
from gene_model import send_mutation, replace_characters, get_reference_from_pmid_by_metapub, pmid_extractor, \
    create_AddLiteratureReferenceJournal_mutation, \
    create_AddLiteratureReferenceAuthors_mutation, fix_author_id, create_author_mutation, create_reference_mutation, \
    create_journal_mutation, get_authors_names, ref_name_from_authors_pmid_and_year, fetch_gene_id_by_gene_name, \
    fetch_gene_info_by_gene_id, populate_omni_gene, create_myGeneInfo_gene, pubmed_extractor, \
    get_gene_id_from_entrez_id, get_omnigene_id_from_entrez_id


def read_edits(file_name:str):
    skip = True;
    edits = []
    with open(file_name) as csvfile:
        readCSV = csv.reader(csvfile, delimiter='\t')
        for row in readCSV:
            if skip:
                skip=False
            else:
                gene_categories = {'gene' : row[0],
                                'description': replace_characters(row[4])}
                edits.append(gene_categories)

    return edits


def send_query(query:str, server:str) -> str:
    url = "http://" + server + ":7474/graphql/"
    headers = {
      'Authorization': 'Basic bmVvNGo6b21uaQ==',
      'Content-Type': 'application/json',
    }
    responseBody = ''
    try:
        response = requests.request("POST", url, headers=headers, json={'query': query})
        if not response.ok:
            response.raise_for_status()
            os.system("say another error")

            sys.exit()

        responseBody: str = response.json()
        # print(responseBody)
        if 'errors' in responseBody:
            print(responseBody)
            os.system("say another error")
            sys.exit()
    except requests.exceptions.RequestException as err:
        print('error in request')
        print(err)
        print(responseBody)
        os.system("say another error")
        sys.exit()
    return responseBody

def get_id_old_id_and_field(gene_name, server):
    id, old_id, field = None,None,None
    query1 = '''{
  OmniGene(name:"'''
    query = query1 + gene_name + '''"){
    id
    geneDescription{
      id
      field
    }
  }
}'''
    response = send_query(query, server)
    if len(response['data']['OmniGene']) > 0:
        result = response['data']['OmniGene'][0]
        id = result['id']
        geneDescription = result['geneDescription']
        old_id = geneDescription['id']
        field = geneDescription['field']
    return id, old_id, field

# {
#   Journal(name:"Science"){
#     id
#   }
# }
def get_journal_id(journal: str,server):
    id = None
    # query1 = '''
    # {
    #   Journal(name:"'''
    # query = query1 + journal + '''"){
    #     id
    #    }
    # }'''
    query = f'{{ Journal(name:"{journal}"){{id}} }}'

    response = send_query(query, server)
    if len(response['data']['Journal'])>0:
        result = response['data']['Journal'][0]
        id = result['id']
    return id

def get_author_id(surname:str,first_initial:str,server):
    id = None
    query = f'{{ Author(surname:"{surname}",first_initial:"{first_initial}"){{id}} }}'

    response = send_query(query, server)
    if len(response['data']['Author'])>0:
        result = response['data']['Author'][0]
        id = result['id']
    return id

def get_reference_by_pmid(pmid:str,server)->str:
    id = None
    ref_id = 'ref_' + pmid
    query = f'{{ LiteratureReference ( id:"{ref_id}" ) {{ id, PMID }} }}'
    response = send_query(query, server)
    if len(response['data']['LiteratureReference'])>0:
        result = response['data']['LiteratureReference'][0]
        id = result['id']
    return id


# addEditableStatementReferences(id: ID!references: [ID!]!): String
# Adds References to EditableStatement entity

def write_references(es_id:str, description:str,server:str )->str:
    s: str = ''
    references = []
    pmids = pmid_extractor(description)
    for pmid in pmids:
        reference = get_reference_from_pmid_by_metapub(pmid)
        if reference != None:
            references.append(reference)
    if len(references)>0:
        reference_string = '['
        for r in references:
            pubmed = r['pmid']
            ref_id = get_reference_by_pmid(pubmed, server)
            if ref_id == None:
                ref_id = 'ref_' + pubmed
                s += create_reference_mutation(ref_id, r)
                journal = r['journal']
                journal_id = get_journal_id(journal, server)
                if journal_id == None:
                    journal_id = 'journal_' + fix_author_id(journal)
                    s += create_journal_mutation(journal, journal_id)
                s += create_AddLiteratureReferenceJournal_mutation(ref_id, journal_id)
                authors = []
                for author in r['authors']:
                    first, surname = get_authors_names(author)
                    author_id = get_author_id(surname, first, server)
                    if author_id == None:
                        author_id = 'author_' + surname + '_' + first
                        author_id = fix_author_id(author_id)
                        s += create_author_mutation(author_id, surname, first)
                    authors.append(author_id)
                s += create_AddLiteratureReferenceAuthors_mutation(ref_id, authors)
            reference_string += '\\"' + ref_id + '\\",'
        reference_string += ']'
        s += f'addEditableStatementReferences(id:\\"{es_id}\\", references:{reference_string}),'
    return s


def createEditableStatement(statement:str, field:str, editor:str,server) -> (str,str):
    now = datetime.datetime.now()
    edit_date:str = now.strftime("%Y-%m-%d-%H-%M-%S-%f")
    id:str = 'es_' + now.strftime("%Y%m%d%H%M%S%f")
    s = f'''{id} : createEditableStatement(deleted: false, edit_date: \\"{edit_date}\\", editor: \\"{editor}\\",field: \\"{field}\\", id: \\"{id}\\",statement: \\"{statement}\\"),'''
    s += write_references(id,statement,server)
    return s, id


# deleteOmniGeneGeneDescription(geneDescription: [ID!]!id: ID!): String
# Deletes GeneDescription from OmniGene entity

# addOmniGeneGeneDescription(geneDescription: [ID!]!id: ID!): String
# Adds GeneDescription to OmniGene entity

def write_new_gene_description(gene_id, old_description_id:str, field:str, statement:str, editor, server:str):
    s: str = '{"query":"mutation {'

    s += f'''deleteOmniGeneGeneDescription(geneDescription:[\\"{old_description_id}\\"], id:\\"{gene_id}\\"),'''
    (m, id1) = createEditableStatement(statement, field, editor,server)
    s += m
    s += f'addOmniGeneGeneDescription(geneDescription:[\\"{id1}\\"], id:\\"{gene_id}\\" ),'
    s += '}"}'
    print(s)
    responseBody = send_mutation(s, server)
    print(responseBody)

def write_uniprot_references(uniprot_id:str, description:str,server:str)->str:
    s: str = ''
    pub_med_ids = pubmed_extractor(description)
    references: list = []
    for pmid in pub_med_ids:
        ref = get_reference_from_pmid_by_metapub(pmid)
        references.append(ref)

    if len(references)>0:
        reference_string = '['
        for r in references:
            pubmed = r['pmid']
            ref_id = get_reference_by_pmid(pubmed,server)
            if ref_id==None:
                ref_id = 'ref_' + pubmed
                s += create_reference_mutation(ref_id, r)
                journal = r['journal']
                journal_id = get_journal_id(journal,server)
                if journal_id == None:
                    journal_id = 'journal_' + fix_author_id(journal)
                    s += create_journal_mutation(journal, journal_id)
                s += create_AddLiteratureReferenceJournal_mutation(ref_id, journal_id)
                authors = []
                for author in r['authors']:
                    first, surname = get_authors_names(author)
                    author_id = get_author_id(surname,first, server)
                    if author_id == None:
                        author_id = 'author_' + surname + '_' + first
                        author_id = fix_author_id(author_id)
                        s += create_author_mutation(author_id, surname, first)
                    authors.append(author_id)
                s += create_AddLiteratureReferenceAuthors_mutation(ref_id, authors)

            reference_string += '\\"' + ref_id + '\\",'
        reference_string += ']'
        s += f'addUniprot_EntryReferences(id:\\"{uniprot_id}\\", references:{reference_string}),'
    return s

def create_uniprot_entry(omni_gene:dict,server:str):
    if 'sp_info' in omni_gene:
        sp_info = omni_gene['sp_info']
        id: str = sp_info['id']
        accessionNumber: str = sp_info['acc_num']
        statement: str = replace_characters(sp_info['function'])
        name: str = sp_info['name']
        uniprot_id: str = sp_info['uniprot_id']
        mutation_payload: str = '{"query":"mutation {'
        s = f'{id}: createUniprot_Entry(accessionNumber: \\"{accessionNumber}\\", id: \\"{id}\\", name: \\"{name}\\", statement: \\"{statement}\\", uniprot_id:\\"{uniprot_id}\\"),'
        mutation_payload += s
# addUniprot_EntryGene(gene: [ID!]!id: ID!): String
# Adds Gene to Uniprot_Entry entity
        gene_id = get_gene_id_from_entrez_id(omni_gene['entrez_gene_id'])
        s = f'addUniprot_EntryGene(gene:[\\"{gene_id}\\"], id:\\"{id}\\" ),'
        mutation_payload += s
        s = write_uniprot_references(id,statement,server)
        mutation_payload += s
        mutation_payload += '}"}'
        send_mutation(mutation_payload, server)


def create_omniGene(omni_gene:dict, gene_dict:dict, gene_description:str, editor, server:str)->None:
    id = get_omnigene_id_from_entrez_id(omni_gene['entrez_gene_id'])
    gene: str = omni_gene['symbol']
    panel_name = omni_gene['panel_name']
    mutation_payload: str = '{"query":"mutation {'
    s = f'{id}: createOmniGene(id: \\"{id}\\", name: \\"{gene}\\", panelName:\\"{panel_name}\\" ),'
    # create geneDescription EditableStatement
    field1: str = 'geneDescription_' + id
    if gene_description==None:
        gene_description = '(Insert Gene Description)'
    statement1: str = gene_description
    (m, id1) = createEditableStatement(statement1,field1,editor,server)
    s += m
#     addOmniGeneGeneDescription(geneDescription: [ID!]!id: ID!): String
#     Adds GeneDescription to OmniGene entity
    s += f'addOmniGeneGeneDescription(geneDescription:[\\"{id1}\\"], id:\\"{id}\\" ),'
    if 'category' in omni_gene:
        # create OncogenicCategory EditableStatement
        field2: str = 'OncogenicCategory_' + id
        statement2: str = '(Draft) ' + omni_gene['category']
        (m, id2) = createEditableStatement(statement2,field2,editor,server)
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
    (m, id3) = createEditableStatement(statement3, field3, editor,server)
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

    # addOmniGeneUniprot_entry(id: ID!uniprot_entry: [ID!]!): String
    # Adds Uniprot_entry to OmniGene entity
    if 'sp_info' in omni_gene:
        uniprot_id:str = omni_gene['sp_info']['id']
        s += f'addOmniGeneUniprot_entry(id:\\"{id}\\", uniprot_entry:[\\"{uniprot_id}\\"] ),'
    mutation_payload: str = '{"query":"mutation {'
    mutation_payload += s + '}"}'
    # print(mutation_payload)
    send_mutation(mutation_payload, server)



def create_omni_gene(gene_name, gene_description, editor,server):
    omni_gene: dict = {
        'symbol': gene_name,
        'panel_name': gene_name
    }
    entrez_gene_id = fetch_gene_id_by_gene_name(gene_name)
    omni_gene['entrez_gene_id'] = entrez_gene_id
    if entrez_gene_id is None:
        print("no entrz gene id for", gene_name)
    else:
        gene_info = fetch_gene_info_by_gene_id(entrez_gene_id)
        populate_omni_gene(gene_info, omni_gene)
        print(omni_gene)
        create_myGeneInfo_gene(omni_gene, server)
        create_uniprot_entry(omni_gene, server)
        create_omniGene(omni_gene,{},gene_description,editor,server)

def main():
    server: str = 'localhost'
    editor = 'Carrie.Hoefer@omniseq.com'
    edits = read_edits('data/GeneDescriptions.tsv')
    for edit in edits:
        print(edit['gene'])
        gene_id, old_description_id,field = get_id_old_id_and_field(edit['gene'], server)
        if gene_id==None:
            print('no gene id for ' + edit['gene'])
            create_omni_gene(edit['gene'],edit['description'],editor, server)
        else:
            write_new_gene_description(gene_id,old_description_id,field,edit['description'], editor,server)
            print(gene_id,old_description_id,field)
        print()




if __name__ == "__main__":
    main()


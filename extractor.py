import csv
import datetime

from graphql_utils import send_query


def get_current_data(server:str)->list:
    query = '''
{
  OmniGene(orderBy: name_asc) {
    name
    geneDescription {
      ...es_fields
    }
    oncogenicCategory {
      ...es_fields
    }
    synonymsString {
      ...es_fields
    }
  }
}
fragment es_fields on EditableStatement {
  statement
  editor {
    name
  }
  edit_date
  references {
    ... on LiteratureReference {
      PMID
      DOI
      title
      journal {
        name
      }
      volume
      first_page
      last_page
      publication_Year
      abstract
      authors {
        surname
        first_initial
      }
    }
  }
}
    '''

    response = send_query(query, server)
    return response['data']['OmniGene']

def get_current_data_fields(server:str)->list:
    query = '''
{
  OmniGene(orderBy: name_asc) {
    name
    geneDescription {
      field
    }
    oncogenicCategory {
      field
    }
    synonymsString {
      field
    }
  }
}
    '''
    response = send_query(query, server)
    return response['data']['OmniGene']

def get_history(server:str)->list:
    query = f'''{{EditableStatement{{
    field
    statement
    edit_date
    editor{{
      name
    }}
  }}
}}'''
    response = send_query(query, server)
    return response['data']['EditableStatement']


def write_descriptions(descriptions:list):
    csv_file = "out/descriptions.csv"
    csv_columns = ['gene_name','gene_description','editor','edit_date','references']
    try:
        with open(csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in descriptions:
                writer.writerow(data)
    except IOError:
        print("I/O error")

def write_oncogenic_categories(oncogenic_categories:list):
    csv_file = "out/oncogenic_categories.csv"
    csv_columns = ['gene_name','oncogenicCategory','editor','edit_date']
    try:
        with open(csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in oncogenic_categories:
                writer.writerow(data)
    except IOError:
        print("I/O error")

def write_synonyms(synonyms:list):
    csv_file = "out/synonyms.csv"
    csv_columns = ['gene_name','synonymsString','editor','edit_date']
    try:
        with open(csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in synonyms:
                writer.writerow(data)
    except IOError:
        print("I/O error")

def write_references(references:list):
    csv_file = "out/references.csv"
    csv_columns = ['PMID', 'DOI', 'title','journal_name','volume','first_page','last_page','publication_Year', 'authors_names']
    try:
        with open(csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns, extrasaction='ignore')
            writer.writeheader()
            for data in references:
                writer.writerow(data)
    except IOError:
        print("I/O error")



def extract_current_data(server:str):
    descriptions: list = []
    oncogenic_categories: list = []
    synonyms: list = []
    reference_dict: dict = {}
    response: list = get_current_data(server)
    for entry in response:
        # print(entry)
        gene_description: dict = entry['geneDescription']
        editor = gene_description['editor']['name']
        edit_date = gene_description['edit_date']
        references = gene_description['references']
        pmids: list = []
        for ref in references:
            if 'PMID' in ref:
                pmids.append(ref['PMID'])
                if ref['PMID'] not in reference_dict:
                    ref['journal_name'] = ref['journal']['name']
                    authors_names = ''
                    if 'authors' in ref:
                        if len(ref['authors']) > 2:
                            authors_names = ref['authors'][0]['surname'] + ' et al.'
                        elif len(ref['authors']) == 2:
                            authors_names = ref['authors'][0]['surname'] + ' & ' + ref['authors'][1]['surname']
                        else:
                            authors_names = ref['authors'][0]['surname']
                    ref['authors_names'] = authors_names
                    reference_dict[ref['PMID']] = ref
        description: dict = {'gene_name': entry['name'], 'gene_description': gene_description['statement'],
                             'editor': editor, 'edit_date': edit_date, 'references': pmids}
        # print(description)
        descriptions.append(description)

        category = entry['oncogenicCategory']
        editor = category['editor']['name']
        edit_date = category['edit_date']
        oncogenic_category: dict = {'gene_name': entry['name'], 'oncogenicCategory': category['statement'],
                                    'editor': editor, 'edit_date': edit_date}
        # print(oncogenic_category)
        oncogenic_categories.append(oncogenic_category)

        synonym_entry = entry['synonymsString']
        editor = synonym_entry['editor']['name']
        edit_date = synonym_entry['edit_date']
        synonym: dict = {'gene_name': entry['name'], 'synonymsString': synonym_entry['statement'], 'editor': editor,
                         'edit_date': edit_date}
        # print(synonym)
        synonyms.append(synonym)
    write_descriptions(descriptions)
    write_oncogenic_categories(oncogenic_categories)
    write_synonyms(synonyms)
    write_references(reference_dict.values())


def write_history(history_list:list):
    csv_file = "out/history.csv"
    csv_columns = ['gene_name', 'kind', 'editor_name', 'edit_date', 'statement']
    try:
        with open(csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns, extrasaction='ignore')
            writer.writeheader()
            for data in history_list:
                writer.writerow(data)
    except IOError:
        print("I/O error")

# 'geneDescription_omnigene_7249'

def convert_field(field:str):
    f = 'GeneDescription'
    if 'SynonymsString' in field:
        f = 'SynonymsString'
    elif 'OncogenicCategory' in field:
        f = 'OncogenicCategory'
    return f

def extract_history(server:str):
    response = get_current_data_fields(server)
    history_list:list = []
    gene_dict = {}

    for entry in response:
        gene_description: dict = entry['geneDescription']
        oncogenicCategory: dict = entry['oncogenicCategory']
        synonymsString: dict = entry['synonymsString']
        field = gene_description['field']
        gene_dict[field] = entry['name']
        field = oncogenicCategory['field']
        gene_dict[field] = entry['name']
        field = synonymsString['field']
        gene_dict[field] = entry['name']

    history_response = get_history(server)
    for h_item in history_response:
        field = h_item['field']
        if field in gene_dict:
            gene = gene_dict[field]
            kind = convert_field(field)
            editor_name = h_item['editor']['name']
            history = {'gene_name': gene, 'kind':kind, 'statement':h_item['statement'], 'editor_name':editor_name, 'edit_date': h_item['edit_date']}
            history_list.append(history)
    write_history(history_list)


def main():
    server:str = 'localhost'
    # server: str = '165.227.89.140'

    extract_current_data(server)
    extract_history(server)



if __name__ == "__main__":
    main()

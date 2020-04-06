import sys
import os
import json
import csv
import requests
from neo4j import GraphDatabase
import typing
import datetime
import time
os.environ['NCBI_API_KEY'] = 'cde5c1a63fa16711994bfe74b858747cbb08'
from metapub import PubMedFetcher
import re

JAX_PATH: str = '/Users/mglynias/jax_march_2020/api-export/'

def send_mutation(mutation_payload:str, server:str) -> str:
    url = "http://" + server + ":7474/graphql/"
    headers = {
      'Authorization': 'Basic bmVvNGo6b21uaQ==',
      'Content-Type': 'application/json',
    }
    responseBody = ''
    try:
        response = requests.request("POST", url, headers=headers, data = mutation_payload)
        if not response.ok:
            response.raise_for_status()
            print(mutation_payload)
            os.system("say another error")

            sys.exit()

        responseBody: str = response.json()
        # print(responseBody)
        if 'errors' in responseBody:
            print(mutation_payload)
            print(responseBody)
            os.system("say another error")
            sys.exit()
    except requests.exceptions.RequestException as err:
        print('error in request')
        print(err)
        print(mutation_payload)
        print(responseBody)
        os.system("say another error")
        sys.exit()
    except UnicodeEncodeError as err:
        print('UnicodeEncodeError')
        print(err)
        print(mutation_payload)
        print(responseBody)
        os.system("say another error")
        sys.exit()
    # print(responseBody)
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
        a_string = a_string.replace(u"\u0391", 'A')
        a_string = a_string.replace('β', 'beta')
        a_string = a_string.replace('γ', 'gamma')
        a_string = a_string.replace('δ', 'delta')
        a_string = a_string.replace('ε', 'epsilon')
        a_string = a_string.replace('ζ', 'zeta')
        a_string = a_string.replace('η', 'eta')
        a_string = a_string.replace('θ', 'theta')
        a_string = a_string.replace('ι', 'iota')
        a_string = a_string.replace('ɩ', 'iota')
        a_string = a_string.replace('κ', 'kappa')
        a_string = a_string.replace('λ', 'lamda')
        a_string = a_string.replace('μ', 'mu')
        a_string = a_string.replace('ν', 'nu')
        a_string = a_string.replace('π', 'pi')
        a_string = a_string.replace('ρ', 'rho')
        a_string = a_string.replace('σ', 'sigma')
        a_string = a_string.replace('χ', 'chi')

        a_string = a_string.replace('ω', 'omega')
        a_string = a_string.replace(u"\u0394", 'delta')

        a_string = a_string.replace(u"\u03c5", 'upsilon')
        a_string = a_string.replace(u"\u03a5", 'Upsilon')
        a_string = a_string.replace('Ψ', 'Psi')
        a_string = a_string.replace('Ω', 'Omega')

        a_string = a_string.replace(u"\u025b", 'e')
        a_string = a_string.replace(u"\u0190", 'e')
        a_string = a_string.replace(u"\u223c", '~')
        a_string = a_string.replace(u"\u301c", '~')
        a_string = a_string.replace("č", 'c')
        a_string = a_string.replace("ć", 'c')
        a_string = a_string.replace("ş", "s")
        a_string = a_string.replace("ś", "s")
        a_string = a_string.replace("š", "s")
        a_string = a_string.replace("Š", "S")
        a_string = a_string.replace("ő", "o")
        a_string = a_string.replace("õ", "o")
        a_string = a_string.replace("ń", "n")
        a_string = a_string.replace("ň", "n")
        a_string = a_string.replace("æ", "ae")
        a_string = a_string.replace("ě", "e")
        a_string = a_string.replace("ė", "e")
        a_string = a_string.replace("ę", "e")
        a_string = a_string.replace("Ş", "S")
        a_string = a_string.replace("ů", "u")
        a_string = a_string.replace("ř", "r")
        a_string = a_string.replace("ﬁ", "fi")
        a_string = a_string.replace("ż", "z")
        a_string = a_string.replace("ź", "z")
        a_string = a_string.replace("ğ", "g")
        a_string = a_string.replace("ť", "t")
        a_string = a_string.replace("ž", "z")
        a_string = a_string.replace("ą", "a")

        a_string = a_string.replace("’", "")
        a_string = a_string.replace('"', '')
        a_string = a_string.replace('\\', ' ')
        a_string = a_string.replace(u"\u2216", ' ')

        a_string = a_string.replace(u"\u201c", '')
        a_string = a_string.replace(u"\u201d", '')
        a_string = a_string.replace(u"\u2018", '')
        a_string = a_string.replace(u"\u2019", '')
        a_string = a_string.replace(u"\u05f3", '')
        a_string = a_string.replace(u"\u2032", '_')
        a_string = a_string.replace(u"\u2020", '*')
        a_string = a_string.replace(u"\u0142", '')
        a_string = a_string.replace(u"\u202f", ' ')
        a_string = a_string.replace(u"\u200a", ' ')
        a_string = a_string.replace(u"\u2002", ' ')
        a_string = a_string.replace('→', '->')
        a_string = a_string.replace(u"\u2012", '-')
        a_string = a_string.replace(u"\u207b", '-')
        a_string = a_string.replace(u"\uff0c", ', ')

        a_string = a_string.replace(u"\u207a", '+')
        a_string = a_string.replace(u"\u2011", '-')
        a_string = a_string.replace(u"\u2013", '-')
        a_string = a_string.replace(u"\u2014", '-')
        a_string = a_string.replace(u"\u2044", '/')
        a_string = a_string.replace(u"\u2122", 'TM')
        a_string = a_string.replace(u"\u2005", ' ')
        a_string = a_string.replace(u"\u2009", ' ')
        a_string = a_string.replace(u"\u0131", 'i')
        a_string = a_string.replace(u"\u2081", '1')
        a_string = a_string.replace(u"\u2082", '2')
        a_string = a_string.replace(u"\u2265", '>=')
        a_string = a_string.replace(u"\u2264", '<=')
        a_string = a_string.replace(u"\u2264", '<=')
        a_string = a_string.replace(u"\u226b", ' >>')
        a_string = a_string.replace(u"\u2248", ' =')
        a_string = a_string.replace('\t',' ')
        a_string = a_string.replace('\r','')
        a_string = a_string.replace('\n','')
        a_string = a_string.replace('⁸⁸','')
        a_string = a_string.replace('⁹⁰','')
        a_string = a_string.replace('Ⅱ','II')
        a_string = a_string.replace('Ⅰ','I')

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
                        # 'references': references,
                        'entrezId' :entrezId,
                        'synonyms': synonyms,
                        'canonicalTranscript': canonicalTranscript,
                        'chromosome': chromosome}
    return gene_dict


def get_reference_from_pmid_by_metapub(pmid:str)->dict:
    fetch = PubMedFetcher()
    reference = None
    try:
        time.sleep(0.1)
        article = fetch.article_by_pmid(pmid)
        reference = {'journal':article.journal,
                     'authors': article.authors,
                     'issue':article.issue,
                     'first_page':article.first_page,
                     'last_page': article.last_page,
                     'volume':article.volume,
                     'year': str(article.year),
                     'abstract': replace_characters(article.abstract),
                     'title': replace_characters(article.title),
                     'doi': article.doi,
                     'pmid': article.pmid
                     }
    except:
        print('*** Bad PMID:',pmid)

    return reference

def pmid_extractor(text:str)->list:
    pattern = r'PMID:\s+\d{8}'
    matches = re.findall(pattern,text)
    pmids = []
    for match in matches:
        if match not in pmids:
            pmids.append(match)
    return pmids

def pubmed_extractor(text:str)->list:
    pattern = r'PubMed:\d{8}'
    matches = re.findall(pattern,text)
    pmids = []
    for match in matches:
        match = match[7:]
        if match not in pmids:
            pmids.append(match)
    return pmids

# createLiteratureReference(
# DOI: String
# PMID: String
# abstract: String
# first_page: String!
# id: ID!
# last_page: String!
# publication_Year: String!
# shortReference: String!
# title: String!
# volume: String!): String
# Creates a LiteratureReference entity


def create_AddLiteratureReferenceJournal_mutation(ref_id, journal_id):
    id = ref_id + '_' + journal_id
    s = f'{id}: addLiteratureReferenceJournal(id:\\"{ref_id}\\", journal:\\"{journal_id}\\"),'
    return s


# addLiteratureReferenceAuthors(authors: [ID!]!id: ID!): String
# Adds Authors to LiteratureReference entity
def create_AddLiteratureReferenceAuthors_mutation(ref_id, authors):
    id = 'author_' +ref_id
    author_string = '['
    for a in authors:
        if len(author_string)>1:
            author_string += ","
        author_string += '\\"' + a + '\\"'
    author_string += ']'
    s = f'{id}: addLiteratureReferenceAuthors(id:\\"{ref_id}\\", authors:{author_string}),'

    return s

# createLiteratureReference(
# DOI: String
# PMID: String
# abstract: String
# first_page: String!
# id: ID!
# last_page: String!
# publication_Year: String!
# shortReference: String!
# title: String!
# volume: String!): String
# Creates a LiteratureReference entity
def create_reference_mutation(ref_id, ref):
    ref_name = ref_name_from_authors_pmid_and_year(ref['authors'], ref['pmid'], ref['year'])
    s = f'''{ref_id}: createLiteratureReference(id: \\"{ref_id}\\", abstract: \\"{ref['abstract']}\\", shortReference: \\"{ref_name}\\", title: \\"{ref['title']}\\", volume: \\"{ref['volume']}\\", first_page: \\"{ref['first_page']}\\", last_page: \\"{ref['last_page']}\\", publication_Year: \\"{ref['year']}\\", DOI: \\"{ref['doi']}\\", PMID: \\"{ref['pmid']},\\"),'''
    return s

# createAuthor(
# first_initial: String!
# id: ID!
# surname: String!): String
# Creates a Author entity
def create_author_mutation(id,surname,first):
    s = f'''{id}: createAuthor(first_initial: \\"{first}\\" , id: \\"{id}\\",surname: \\"{surname}\\"),'''
    return s


def create_journal_mutation(journal, journal_id):
    s = f'''{journal_id}: createJournal(id: \\"{journal_id}\\",name: \\"{journal}\\"),'''
    return s

def get_authors_names(author):
    l = author.split()
    surname = replace_characters(l[0])
    first = '-'
    if len(l)>1:
        first = replace_characters(l[1])
    return first, surname

def ref_name_from_authors_pmid_and_year(authors, pmid, year):
    s = ''
    if len(authors)>0:
        first, surname = get_authors_names(authors[0])
        if len(authors) == 1:
            s += surname + ' ' + year
        elif len(authors) == 2:
            first2, surname2 = get_authors_names(authors[1])
            s += surname + '& '+ surname2 + ' ' + year
        else:
            s += surname + ' et al. ' + year
    else:
        s += 'no_authors ' + year
    s += ' (PMID:' + pmid + ')'
    return s


def get_gene_id_from_jax_id(jax_id:str)->str:
    return 'jax_gene_' +jax_id



def fix_author_id(id:str)->str:
    id = id.lower()
    id = id.replace(" ", "")
    id = id.replace(",", "")
    id = id.replace("(", "")
    id = id.replace(")", "")
    id = id.replace("<sup>®<_sup>","")
    id = id.replace("<", "")
    id = id.replace(">", "")
    id = id.replace("®", "_")
    id = id.replace("-", "_")
    id = id.replace("'", "_")
    id = id.replace(".", "_")
    id = id.replace("/", "_")
    id = id.replace("á", "a")
    id = id.replace("à", "a")
    id = id.replace("ä", "a")
    id = id.replace("å", "a")
    id = id.replace("ã", "a")
    id = id.replace("â", "a")


    id = id.replace("æ", "ae")

    id = id.replace("ç", "c")
    id = id.replace("č", "c")

    id = id.replace("é", "e")
    id = id.replace("è", "e")
    id = id.replace("ë", "e")
    id = id.replace("ě", "e")
    id = id.replace("ê", "e")

    id = id.replace("í", "i")
    id = id.replace("ì", "i")
    id = id.replace("î", "i")
    id = id.replace("ï", "i")
    id = id.replace("ñ", "n")
    id = id.replace("ń", "n")
    id = id.replace("ö", "o")
    id = id.replace("ó", "o")
    id = id.replace("ò", "o")
    id = id.replace("ő", "o")
    id = id.replace("ô", "o")
    id = id.replace("ø", "o")
    id = id.replace("ş", "s")
    id = id.replace("ß", "s")
    id = id.replace("ü", "u")
    id = id.replace("ú", "u")
    id = id.replace("ÿ", "y")
    id = id.replace("ý", "y")

    return id

def write_jax_references(gene_id:str, description:str,reference_dict:dict, journal_dict:dict, author_dict:dict )->str:
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
            if not pubmed in reference_dict:
                ref_id = 'ref_' + pubmed
                reference_dict[pubmed] = ref_id
                s += create_reference_mutation(ref_id, r)
                journal = r['journal']
                if not journal in journal_dict:
                    journal_id = 'journal_' + fix_author_id(journal)
                    journal_dict[journal] = journal_id
                    s += create_journal_mutation(journal, journal_id)
                else:
                    journal_id = journal_dict[journal]
                s += create_AddLiteratureReferenceJournal_mutation(ref_id, journal_id)
                authors = []
                for author in r['authors']:
                    first, surname = get_authors_names(author)
                    id = 'author_' + surname + '_' + first
                    id = fix_author_id(id)
                    if not id in author_dict:
                        author_dict[id] = id
                        s += create_author_mutation(id, surname, first)
                        authors.append(id)
                s += create_AddLiteratureReferenceAuthors_mutation(ref_id, authors)
            else:
                ref_id = reference_dict[pubmed]
            reference_string += '\\"' + ref_id + '\\",'
        reference_string += ']'
        s += f'addJaxGeneReferences(id:\\"{gene_id}\\", references:{reference_string}),'
    return s

# createJaxGene(
# canonicalTranscript: [String]!
# chromosome: String!
# entrezId: String!
# id: ID!
# jaxId: String!
# name: String!
# statement: String!
# synonyms: [String]!): String

def write_one_jax_gene(gene: dict,gene_dict:dict,reference_dict:dict,journal_dict: dict,author_dict: dict)->str:
    id: str = get_gene_id_from_jax_id(str(gene['ckb_id']))
    gene_dict[gene['gene']] = id
    synonyms: str = '['
    for syn in gene['synonyms']:
        synonyms += f'\\"{syn}\\",'
    synonyms += ']'

    s = f'''{id} : createJaxGene(canonicalTranscript: \\"{gene['canonicalTranscript']}\\", chromosome: \\"{gene['chromosome']}\\",entrezId: \\"{gene['entrezId']}\\", id: \\"{id}\\", jaxId: \\"{gene['ckb_id']}\\", name: \\"{gene['gene']}\\", statement: \\"{gene['description']}\\",synonyms: {synonyms}),'''
    s += write_jax_references(id,gene['description'],reference_dict,journal_dict,author_dict)
    return s


def write_jax_genes(server:str,gene_dict:dict, reference_dict: dict,journal_dict: dict,author_dict: dict) -> None:
    # gene_categories = read_oncgenes_tumor_suppressors('data/OCP_oncogenes_tumorsuppressors.csv')
    json_files = get_list_of_files(JAX_PATH + 'genes')
    print('total=',len(json_files))
    count = 0
    last_good = -1
    for json_file in json_files:
        gene: dict = read_one_gene_jax_json(json_file)
        if gene is not None:
            if count > last_good:
                mutation_payload: str = '{"query":"mutation {'
                mutation_payload += write_one_jax_gene(gene,gene_dict,reference_dict,journal_dict,author_dict)
                mutation_payload += '}"}'
                print(count,gene['gene'])
                send_mutation(mutation_payload, server)
            count += 1



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


def create_mygene_reference(gene_id):
  reference = {}
  reference['type'] = 'InternetReference'
  reference['url'] = "http://mygene.info/v3/gene/" + gene_id
  reference['accessed_date'] = datetime.datetime.now()
  return reference


def fetch_uniprot_by_acc_num(id):
  requestURL = "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession=" + id
  r = requests.get(requestURL, headers={ "Accept" : "application/json"})
  if not r.ok:
    r.raise_for_status()
    sys.exit()
  responseBody = r.json()
  # pprint.pprint(responseBody)
  return responseBody[0]


def get_sp_info(uniprot_accession:str)->dict:
    sp_info = {'function': '', 'name':''}

    sp: dict = fetch_uniprot_by_acc_num(uniprot_accession)
    sp_info['acc_num'] = sp['accession']
    sp_info['id'] = 'uniprot_' + sp['accession']
    sp_info['uniprot_id'] = sp['id']
    if 'recommendedName' in sp['protein']:
        sp_info['name'] = sp['protein']['recommendedName']['fullName']['value']
    elif 'submittedName' in sp['protein']:
        sp_info['name'] = sp['protein']['submittedName'][0]['fullName']['value']

    for comment in sp['comments']:
        if comment['type'] == 'FUNCTION':
            if len(comment['text'])>0 and 'value' in comment['text'][0]:
                sp_info['function'] = comment['text'][0]['value']
            break

    return sp_info

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

def createEditableStatement(statement:str, field:str, editor:str) -> (str,str):
    now = datetime.datetime.now()
    edit_date:str = now.strftime("%Y-%m-%d-%H-%M-%S-%f")
    id:str = 'es_' + now.strftime("%Y%m%d%H%M%S%f")
    s = f'''{id} : createEditableStatement(deleted: false, edit_date: \\"{edit_date}\\", editor: \\"{editor}\\",field: \\"{field}\\", id: \\"{id}\\",statement: \\"{statement}\\"),'''
    return s, id

# addUniprot_EntryReferences(id: ID!references: [ID!]!): String
# Adds References to Uniprot_Entry entity
def write_uniprot_references(uniprot_id:str, description:str,reference_dict:dict, journal_dict:dict, author_dict:dict )->str:
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
            if not pubmed in reference_dict:
                ref_id = 'ref_' + pubmed
                reference_dict[pubmed] = ref_id
                s += create_reference_mutation(ref_id, r)
                journal = r['journal']
                if not journal in journal_dict:
                    journal_id = 'journal_' + fix_author_id(journal)
                    journal_dict[journal] = journal_id
                    s += create_journal_mutation(journal, journal_id)
                else:
                    journal_id = journal_dict[journal]
                s += create_AddLiteratureReferenceJournal_mutation(ref_id, journal_id)
                authors = []
                for author in r['authors']:
                    first, surname = get_authors_names(author)
                    id = 'author_' + surname + '_' + first
                    id = fix_author_id(id)
                    if not id in author_dict:
                        author_dict[id] = id
                        s += create_author_mutation(id, surname, first)
                    authors.append(id)
                s += create_AddLiteratureReferenceAuthors_mutation(ref_id, authors)
            else:
                ref_id = reference_dict[pubmed]
            reference_string += '\\"' + ref_id + '\\",'
        reference_string += ']'
        s += f'addUniprot_EntryReferences(id:\\"{uniprot_id}\\", references:{reference_string}),'
    return s

# createUniprot_Entry(
# accessionNumber: String!
# id: ID!
# name: String!
# statement: String!
# uniprot_id: String!): String
# Creates a Uniprot_Entry entity
def create_uniprot_entry(omni_gene:dict,server:str,reference_dict,journal_dict,author_dict):
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
        s = write_uniprot_references(id,statement,reference_dict,journal_dict,author_dict)
        mutation_payload += s
        mutation_payload += '}"}'
        send_mutation(mutation_payload, server)

# createOmniGene(
# id: ID!
# name: String!
# panelName: String!): String
# Creates a OmniGene entity
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
    (m, id1) = createEditableStatement(statement1,field1,editor)
    s += m
#     addOmniGeneGeneDescription(geneDescription: [ID!]!id: ID!): String
#     Adds GeneDescription to OmniGene entity
    s += f'addOmniGeneGeneDescription(geneDescription:[\\"{id1}\\"], id:\\"{id}\\" ),'
    # create OncogenicCategory EditableStatement
    field2: str = 'OncogenicCategory_' + id
    statement2: str = '(Draft) ' + omni_gene['category']
    (m, id2) = createEditableStatement(statement2,field2,editor)
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
    (m, id3) = createEditableStatement(statement3, field3, editor)
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


def write_omni_genes(omini_genes:list, gene_dict:dict, server:str,reference_dict,journal_dict,author_dict):
    hgnc_gene_name_dict = create_hgnc_gene_name_dict()
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
                populate_omni_gene(gene_info, omni_gene)
                print(gene_name)
                create_myGeneInfo_gene(omni_gene,server)
                create_uniprot_entry(omni_gene,server,reference_dict,journal_dict,author_dict)
                create_omniGene(omni_gene,gene_dict,None,'auto',server)
            else:
                print('no gene_info for',gene_name)


def populate_omni_gene(gene_info, omni_gene):
    genomic_pos = gene_info['genomic_pos_hg19']
    if isinstance(genomic_pos, list):
        genomic_pos = genomic_pos[0]
    omni_gene['chrom'] = genomic_pos['chr']
    omni_gene['strand'] = 'FORWARD'
    if genomic_pos['strand'] == -1:
        omni_gene['strand'] = 'REVERSE'
    omni_gene['start'] = genomic_pos['start']
    omni_gene['end'] = genomic_pos['end']
    omni_gene['reference'] = create_mygene_reference(omni_gene['entrez_gene_id'])
    if 'summary' in gene_info:
        omni_gene['summary'] = gene_info['summary']
    if 'alias' in gene_info:
        omni_gene['synonyms'] = gene_info['alias']
    omni_gene['synonyms'] = []
    if 'uniprot' in gene_info:
        if 'Swiss-Prot' in gene_info['uniprot']:
            if isinstance(gene_info['uniprot']['Swiss-Prot'], list):
                uniprot_accession = gene_info['uniprot']['Swiss-Prot'][0]
            else:
                uniprot_accession = gene_info['uniprot']['Swiss-Prot']
            omni_gene['sp_info'] = get_sp_info(uniprot_accession)


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


def main():
    server:str = 'localhost'
    # server: str = '165.227.89.140'
    erase_neo4j(server)

    gene_dict: dict = {}
    reference_dict: dict = {}
    journal_dict: dict = {}
    author_dict: dict = {}

    write_jax_genes(server,gene_dict,reference_dict,journal_dict,author_dict)
    omini_genes:list = read_omni_genes('data/tblOS_GLOBAL_GLOBAL_Ref_AllPanelsGenes.csv')
    write_omni_genes(omini_genes,gene_dict,server,reference_dict,journal_dict,author_dict)



if __name__ == "__main__":
    main()
#
# Define function to download gzipped file from url, unzip it, and load it into pandas DF
def download_and_load_data(url, output_filename, sep="\t", column_filter=None):
    # Step 1: Download the file using wget
    print(f"Downloading {output_filename} from {url}...")
    subprocess.run(["wget", "-O", output_filename + ".gz", url], check=True)

    # Step 2: Gunzip the file
    print(f"Unzipping {output_filename}.gz...")
    with gzip.open(output_filename + ".gz", "rb") as gz_file:
        with open(output_filename, "wb") as out_file:
            out_file.write(gz_file.read())

    # Step 3: Load the data into a pandas DataFrame
    print(f"Loading {output_filename} into a pandas DataFrame...")
    df = pd.read_csv(output_filename, sep=sep, index_col=0)

    # Optional: Filter columns based on the keyword
    if column_filter:
        print(f"Filtering columns with keyword '{column_filter}'...")
        filtered_columns = [col for col in df.columns if column_filter in col]
        df = df[filtered_columns]

    return df

# Load data for 5xFAD mouse model (8mo only)
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168137&format=file&file=GSE168137%5FcountList%2Etxt%2Egz"
output_filename = "GSE168137_countList.txt"
column_keyword = "cortex_4mon"
countlist_5xFAD = download_and_load_data(url, output_filename, column_filter=column_keyword)

# Rename columns
new_column_names = ["5xFAD_cortex_4mo_Female_1", "5xFAD_cortex_4mo_Female_2","5xFAD_cortex_4mo_Female_3", "5xFAD_cortex_4mo_Female_4", "5xFAD_cortex_4mo_Female_5", "5xFAD_cortex_4mon_Male_1", "5xFAD_cortex_4mon_Male_2", "5xFAD_cortex_4mon_Male_3", "5xFAD_cortex_4mon_Male_4", "5xFAD_cortex_4mon_Male_5",  "BL6_cortex_4mon_Female_1", "BL6_cortex_4mon_Female_2", "BL6_cortex_4mon_Female_3",  "BL6_cortex_4mon_Female_4", "BL6_cortex_4mon_Female_5", "BL6_cortex_4mon_Male_1", "BL6_cortex_4mon_Male_2", "BL6_cortex_4mon_Male_3", "BL6_cortex_4mon_Male_4", "BL6_cortex_4mon_Male_5"]
countlist_5xFAD.columns = new_column_names

# Drop ensemble version ID from gene_id's
countlist_5xFAD.index = countlist_5xFAD.index.str.split('.').str[0]

# View first 5 rows of data
countlist_5xFAD.head()

# create MyGeneInfo object
mg = mygene.MyGeneInfo()

# get the ensembl id from index
ensembl_ids = countlist_5xFAD.index.tolist()

# query the gene symbols for the ensemble ids and onvert result to dataframe
gene_info = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='mouse')
gene_df = pd.DataFrame(gene_info)

# remove duplicate ensemble ids and rows where symbol is missing or duplicated
gene_df = gene_df.dropna(subset=['symbol']).drop_duplicates(subset='query')

# map gene symbols back to original dataframe and move gene_name column to front column
countlist_5xFAD['Gene_Name'] = countlist_5xFAD.index.map(gene_df.set_index('query')['symbol'])
cols = ['Gene_Name'] + [col for col in countlist_5xFAD.columns if col != 'Gene_Name']
countlist_5xFAD = countlist_5xFAD[cols]

# view first 5 rows of data
countlist_5xFAD.head()

# rename countlist_5xFAD for use in tutorial 
data = countlist_5xFAD

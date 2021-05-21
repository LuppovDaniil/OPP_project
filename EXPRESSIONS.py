#!/usr/bin/env python
# coding: utf-8

# In[ ]:


class EXPRESSIONS:
    
    def __init__(self, path_rnaseq, path_patient, path_sample):
        self.rna_seq = pd.read_csv(path_rnaseq, sep = "\t")
        self.data_patient = pd.read_csv(path_patient, sep = "\t")
        self.path_sample = pd.read_csv(path_sample, sep = "\t")

    def gene_distribution(self, gene_names):
        '''
        returns genes' distribution histograms
        
        gene_names - array-like object that contains gene names. GeneFinder is used for the genes to be identified.
        '''
        rna_seq = a.rna_seq.copy() # создаём копию, чтобы избежать лишних проблем
        rna_seq.set_index("Entrez_Gene_Id", inplace= True)
        corrected_names = GeneFinder(gene_names)
        plt.figure(figsize=(12, 7))
        for i, j in zip(corrected_names, range(0, len(corrected_names))):
                #print(i, j, rna_seq.loc[i][2::])
                plt.subplot(2, 2, j+1)
                sns.distplot(rna_seq.loc[int(i)][2::], axlabel = rna_seq.loc[int(i)][0])
    
    def clustering(self, Cophenet = True):
        '''
        Calculating the clusterisation of data based on rna sequencing.
        
        Coherent - showing cophenetic correlation of clusterisation
        '''
        
        rna_seq = self.rna_seq.set_index("Hugo_Symbol").copy() #теперь чтобы обратиться к строке нужно использовать имя гена, а не индекс
        #print(rna_seq)
        rna_mad = rna_seq.drop(['Entrez_Gene_Id'], axis = 1).T.mad()
        rna_seq.insert(1, "mad", np.array(rna_mad))
        
        rna_seq_clast1 = rna_seq.copy() #выделим датасет, в котором проведем класетризацию
        rna_seq_clast1 = rna_seq_clast1.sort_values(["mad"], ascending=False) #отберем наиболее альтернативно экспресируемые гены до z-скорирования,
        # тк оно выравнивает все среднеквадратичные отклонения, а метрика mad ей аналогична
        rna_seq_clast1 = rna_seq_clast1[:1500]
        rna_seq_clast1 = rna_seq_clast1.drop("Entrez_Gene_Id", 1)
        rna_seq_clast1 = rna_seq_clast1.drop("mad", 1)
        rna_seq_clast1 = (rna_seq_clast1.T - rna_seq_clast1.T.median())/rna_seq_clast1.T.std()
        rna_seq_clast1 = rna_seq_clast1.dropna(1) #выбросим гены, вернувшие None
        
        
        
        rna_seq_clast_arr1 = np.array(rna_seq_clast1)
        Z1 = sci.linkage(rna_seq_clast_arr1, 'ward')
        
        if Cophenet:
            c, coph_dists = sci.cophenet(Z1, pdist(rna_seq_clast_arr1))
            print(f'Cophenetic correlation is equal to {c}')

        return CLUSTERING(Z1, rna_seq_clast1)
    
    def deconvolution(self, n = 1):
        '''
        returns umap plot of RNA sequence data.
        
        n - the number of clusters. If n != 1, clusters will be shown on the plot.
        '''
        
        reducer = umap.UMAP()
        
        clusters = self.clustering(Cophenet = False) #создаём объект класса CLUSTERING, чтобы использовать таблицу с отпроцессированными экспрессиями
        embedding = reducer.fit_transform(clusters.expressions)
        #print(embedding.shape, embedding)
        
        if n != 1:
            expression_clusterized = clusters.expression_clustering(n) #первый столбец с номерами кластеров
            
            for i in set(expression_clusterized["Cluster"]):
                temp = embedding[expression_clusterized["Cluster"] == i]
                plt.scatter(temp.T[0], temp.T[1])
        else:
            plt.scatter(embedding.T[0], embedding.T[1])

            
def GeneFinder(genes):
    '''
    returns entrezgene ID
    
    genes - array-like object of gene names
    '''
    
    genes_ID = []
    for gene in genes:
        if mg.query(gene, fields='entrezgene')['hits']:
            try:
                genes_ID.append(mg.query(gene, fields='entrezgene')['hits'][0]['entrezgene'])
            except (KeyError):
                warnings.warn('Gene {} not found'.format(gene))
                pass
        else:
            warnings.warn('Gene {} not found'.format(gene))
    return (genes_ID)


#!/usr/bin/env python
# coding: utf-8

# In[ ]:


class CLUSTERING:
    def __init__(self, Z, expressions):
        """
        class for cluterized data
        Z - scoring output
        expressions - pd.dataframe with IDs and expressions
        """
        self.Z = Z
        self.expressions = expressions
        
    def dendogram(self, y = False):
        """
        Build dendogram for hierachial clustering
        param:
        y - cutline
        """
        plt.figure(figsize=(25, 10)) #строим двоичное дерево для нормированных данных
        plt.title('Кластеризация нормированных данных')
        plt.xlabel('sample index')
        plt.ylabel('distance')
        sci.dendrogram(
          self.Z,
          leaf_rotation=9.,  
          leaf_font_size=2., 
        )
        if y:
            plt.axhline(y=y, c='r')
        plt.show()
  
    def expression_clustering(self, n):
        """
        returns expression table with clusters
        param:
        n - number of slusters the data will be clusterized
        """
        clasters = sci.fcluster(self.Z, n, criterion='maxclust')
        expression_clusterized = self.expressions.copy()
        expression_clusterized.insert(0, "Cluster", np.array(clasters))
        return (expression_clusterized)
    
    def sample_clustering(self, n, samples):
        """
        returns sample table with clusters
        param:
        n - number of clusters
        samples - pd.dataframe with samples without service rows
        """
        expression_clusterized = self.expression_clustering(n)
        sample_proc = samples.copy()
        indexes = list(set(sample_proc.index).intersection(set(expression_clusterized.index)))
        sample_proc.insert(0, "Cluster", -1)
        for ind in indexes:
            sample_proc["Cluster"][ind] = expression_clusterized["Cluster"][ind]
        broken = sample_proc[sample_proc["Cluster"] == -1].index
        if len(broken):
            warnings.warn('{} samples were droped'.format(len(broken)))
        sample_proc.drop(broken, inplace = True)
        return (sample_proc)
    
    def patient_clustering(self, n, samples, patients):
        """
        returns patient table with clusters
        param:
        n - number of clusters
        samples - pd.dataframe with samples without service rows
        patients - pd.dataframe with patients without service rows
        """
        sample_proc = self.sample_clustering(n, samples)
        patients_proc = patients.copy()
        patients_proc.insert(0, "Cluster", -1)
        for ind in sample_proc.index:
            patients_proc['Cluster'][sample_proc['#Patient Identifier'][ind]] = sample_proc['Cluster'][ind]
        broken = patients_proc[patients_proc["Cluster"] == -1].index
        if len(broken):
            warnings.warn('{} patient forms were droped'.format(len(broken)))
        patients_proc.drop(broken, inplace = True)
        return patients_proc
    
    def Caplan_Meier(self, n, samples, patients, ci_show = False):
        """
        returns Caplan_Meier diagram with track for each cluster
        param:
        n - number of clusters
        samples - pd.dataframe with samples without service rows
        patients - pd.dataframe with patients without service rows
        """
        patients_proc = self.patient_clustering(n, samples, patients)
        patients_proc.drop(patients_proc[patients_proc['Overall Survival (Months)'] == '[Not Available]'].index, inplace = True)
        events = []
        durations = []
        km = []
        plt.figure(figsize=(25, 10))
        for i in range (n):
            events = list(patients_proc['Overall Survival Status'][patients_proc['Cluster'] == i+1].apply(lambda x: int(x[0])))
            durations = list(patients_proc['Overall Survival (Months)'][patients_proc['Cluster'] == i+1].apply(float))
            kk = KaplanMeierFitter()
            km.append(kk.fit(np.array(durations)/12, events, label='claster {}'.format(i+1)))
        for i in range (n):
            km[i].plot(ci_show=ci_show)
        plt.title('Survival analysis, {} clusters'.format(n), size = 15)
        plt.xlabel('Time, years', size = 15)
        plt.ylabel('Survival rate', size = 15)
        #return (km)
    
    def Heat_Map_Overexpression(self, n, *args):
        """
        returns heatmap with overexpressed genes
        param - number of clusters
        *args - threshold of overexpression for each cluster
        """
        exp = self.expression_clustering(n)
        gene_lists = []
        for i in range (1, n+1):
            exp_i = exp[exp['Cluster'] == i].copy()
            for j in list(exp_i.columns)[1:]:
                if exp_i[j].mean() > exp[j].mean() + thresholds[i-1]:
                    gene_lists.append(j)
        
        plt.figure(figsize=(20, 10))
        plt.title('Upregulation heatmap', size = 15)
        ax = sns.heatmap(exp.sort_values(["Cluster"], ascending=True)[gene_lists], cmap="ocean")
        
    def Heat_Map_Downexpression(self, n, *args):
        """
        returns heatmap with overexpressed genes
        param - number of clusters
        *args - threshold of overexpression for each cluster
        """
        exp = self.expression_clustering(n)
        gene_lists = []
        for i in range (1, n+1):
            exp_i = exp[exp['Cluster'] == i].copy()
            for j in list(exp_i.columns)[1:]:
                if exp_i[j].mean() < exp[j].mean() - thresholds[i-1]:
                    gene_lists.append(j)
        plt.figure(figsize=(20, 10))
        plt.title('Downregulation heatmap', size = 15)
        ax = sns.heatmap(exp.sort_values(["Cluster"], ascending=True)[gene_lists], cmap="ocean")


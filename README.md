# OPP_project_Expression analysis  
Этот проект предоставляет пользователю достаточно исчерпывающий pipeline для анализа экспрессий.  
Для работы с этим проектом пользователь должен предоставить 3 таблицы: таблицу с информацией об образцах, таблицу с инфрмацией о пациентах и таблицу с экспрессиями генов в каждом из образцов.  
## Технологии, используемые в проекте:

- pandas - работа с таблицами
- matplotlib.pyplot - построение графиков
- seaborn - статистика
- numpy - векторы и матрицы
- scipy.cluster.hierarchy - кластеризация
- scipy.spatial.distance - анализ качества кластеризации
- lifelines - анализ выживаемости
- mygene - парсинг генов
- warnings - предупреждения
- umap - деконвалюция методом UMAP 
 
## Функции, реализуемые в этом проекте:
**GeneFinder (genes)**  
returns entrezgene ID
param:
genes - array-like object of gene names  
  
## Классы, реализуемые в этом проекте:  
### - **EXPRESSIONS** 
  #### Инициализация:
  Class for expressions data  
  (path_rnaseq, path_patient, path_sample) - paths for expressions data, patients, samples 
  #### Методы класса:
  - **gene_distribution(self, gene_names)**    
        Returns genes' distribution histograms  
        param:  
        gene_names - array-like object that contains gene names.
  - **clustering(self, Cophenet = True)**  
        Calculating the clusterisation of data based on rna sequencing.  
        param:  
        Coherent - showing cophenetic correlation of clusterisation.  
  - **deconvolution(self, n = 1)**  
        Returns umap plot of RNA sequence data.  
        param:  
        n - the number of clusters. If n != 1, clusters will be shown on the plot.  
### - **CLUSTERING**
  #### Инициализация:
  Сlass for cluterized data  
  (Z, expressions) - scoring output, pd.dataframe with IDs and expressions
  #### Методы класса:
  - **dendogram(self, y = False)**    
        Build dendogram for hierachial clustering  
        param:  
        y - cutline  
  - **expression_clustering(self, n)**  
        Returns expression table with clusters  
        param:  
        n - number of slusters the data will be clusterized  
  - **sample_clustering(self, n, samples)**  
        Returns sample table with clusters  
        param:  
        n - number of clusters  
        samples - pd.dataframe with samples without service rows  
  - **patient_clustering(self, n, samples, patients)**    
        Returns patient table with clusters  
        param:  
        n - number of clusters  
        samples - pd.dataframe with samples without service rows  
        patients - pd.dataframe with patients without service rows  

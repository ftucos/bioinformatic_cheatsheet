# Useful python functions

```python
>>> import os
>>> os.getcwd()
'/Users/tucos'
```

## Data wrangling

```python
>>> import numpy as np
>>> import pandas as pd

>>> df = pd.read_csv(path)
>>> df = pd.read_json(path)
>>> df = pd.read_excel()
>>> ...

>>> df.shape
(1000, 45)

# export
df.to_csv('output.csv',index=False)

# report of number of columns, data type and number of nulls
>>> df.info()

# counts number of nulls for each column
>>> df.isnull().sum()
# fill NAs with 0s
>>> df['selected_column'].fillna(0, inplace=True)
# or fill NAs with the most common value in that column
>>> df['selected_column'].fillna(df['selected_column'].mode, inplace=True)

# computee counts, mean, std and quantiles for eah column
>>> df.describe()

# find the number of unique variables and their count (the same as table() in R)
>>> df['seelcted_column'].value_counts()
# or just get the number of unique values
>>> df['seelcted_column'].nunique()

# drop a column
>>> df.drop(['selected_column'], axis=1, inplace=True)

# get 4 rows with the smallest/largest value for a specific column
>>> df.nsmallest(4, 'selected_column')
>>> df.nlargest(4, 'selected_column')
                
# set a column as index
>>> df.set_index('selected_Column'. inplace=True)

# drop duplicated rows
>>> df.duplicated().sum()
5
>>> df.drop_duplicates(inplace=True)

rbind = pd.concat([df1, df2], axis=0)
cbind = pd.concat([df1,df2], axis=1)

# Convert categorical features into  dummy variables (mask for regression models)
>>> dummy = pd.get_dummies(df.race_column, prefix='race_')

# split dataframes by data types
>>> numerical_df = df.selecte_dtypes(include=[np.number])
>>> categorical_df = df.select_ftypes(include='object')
```

### query()

```python
# filter rows
>>> df.query('45 < `student age` < 55')
>>> df.query('col1 == col2 and col1 != col 3')
# the same as
>>> df.query[(df['col1'] == df['col2']) & (df['col1'] != df['col3'])]

>>> selected_genes = ['ACTB', 'GAPDH', ...]
>>> df.query('genes in @selected_genes')
```

### Summarise data

```python
# summarize method for all the columns
>>> df.groupby('selected_column').agg(np.mean)
# summarize method for a specific column 
>>> df.groupby('selected_column').another_column.describe()
```

### loc() and iloc()

```python
df.loc[:5, ['col1', 'col2']]

# filter rows: [,:] is implicit
df.loc[(df.col1 < 10) & (df.col2 == 'high')]

df.iloc[:4, :5]

```

### apply()

```python
>>> def currentAge(age):
      return (2021 - age)
>>> df['current_age'] = df['birth'].apply(currentAge)
```

### group in bins

```python
# split eavenly into 5 bins
df['height_group'] = pd.qcut(df['height'],q=5)
# define bin sizes
df['height_group'] = pd.cut(df['height'],bins=3, labels=['0-150','151-180','181-300'])



```


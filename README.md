This repository contains impelementation of DNBGDA model[1]

To run the code, you need to prepare two R objects X, and m_len. X is the list of term-document matrix (tm package) and m_len denoces number of documents from each
source at eatch time stamp. You can load the example dataset to get familiar with the required format:

```
load("example_data.RData")
```

After preparing the objects, you can load the required functions by

```
source("DNBGFA_var_functions.R")
```

and then

```
model<- model_setup(X, m_len, 20)
update_model(model, iter = 1000, mc.cores = 10, print_topics = T)

```

or do

```
model<- model_setup(X, m_len, 20)
model<- update_model_interp(model, iter = 10, inter_t = 2, mc.cores = 10)
```

for interpolation task (the second time stamp).

You can also check the file example_run.R to get familiar with the process.

[1] Chien Lu, Jaakko Peltonen, Jyrki Nummenmaa, and Kalervo Järvelin. Probabilistic Dynamic Non-negative Group Factor Model for Multi-source Text Mining. ACM International Conference on Information and Knowledge Management, 2020.

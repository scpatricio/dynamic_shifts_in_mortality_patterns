# Detecting Late-Life Mortality Plateaus  
### Replication Materials for Manuscript III

This folder contains the code accompanying the thesis chapter on **late-life mortality plateaus** (Chapter 3). The chapter investigates whether age-specific mortality stops increasing at the highest ages and develops a framework to **detect and quantify mortality plateaus** using Human Mortality Database (HMD) data.

The material here focuses on:

- estimating age-specific hazards at old ages,
- assessing the deviation from log-linear Gompertz behaviour,
- and identifying ages where mortality levels off into a plateau.

---

## Purpose of This Folder

The goal of this directory is to provide a **reusable toolkit** to:

- read age-specific death counts and exposures from HMD-like files,
- estimate hazards and model them (e.g. Gompertz vs plateau),
- and apply statistical toold to detect wether a **late-life plateau** is present.

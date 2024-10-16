import statsmodels.tsa
import statsmodels.tsa.stattools


noise_space = pd.DataFrame(0, index=facs_2023_mean.columns, columns=["log2mean","CV","mean_pop_CV","T1/2"], dtype=float)
noise_space["log2mean"] = np.log2(facs_2023_mean.mean())
noise_space["CV"] = (facs_2023_mean.std()/facs_2023_mean.mean())
noise_space["mean_pop_CV"] = (facs_2023_cv/100).mean()
acfs = {}
for clone in facs_2023_mean.columns:
    acfs[clone] = statsmodels.tsa.stattools.acf(facs_2023_mean[clone])

for clone, acf in acfs.items():
    for lag, corr in enumerate(acf):
        if corr < 0.5:
            noise_space.loc[clone, "T1/2"] = lag
            break


        
noise_space = noise_space.sort_index()
noise_space

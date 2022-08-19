#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

low_cutoff, hi_cutoff = 0.7, 1.3
higher_than_hi_factor = 1.5

fc_types = {
    "mat_down_pre_undetected": "Mature down, precursor not detected",
    "mat_up_pre_undetected": "Mature up, precursor not detected",
    "mat_undetected_pre_down": "Precursor down, mature not detected",
    "mat_undetected_pre_up": "Precursor up, mature not detected",
    "mat_down_pre_up": "Mature down, precursor up",
    "mat_up_pre_up": "Mature up, precursor up more",
    "mat_down_only": "Only mature down",
    "mat_up_only": "Only mature up",
    "pre_down_only": "Only precursor down",
    "pre_up_only": "Only precursor up",
    "other": "Undetected/Other",
}


def classify_fc(mat_fc, pre_fc):
    if pd.isna(pre_fc) or pre_fc == 0.0:
        if 0.0 < mat_fc <= low_cutoff:
            return fc_types["mat_down_pre_undetected"]
        elif mat_fc >= hi_cutoff:
            return fc_types["mat_up_pre_undetected"]
        return fc_types["other"]
    if pd.isna(mat_fc) or mat_fc == 0.0:
        if 0.0 < pre_fc <= low_cutoff:
            return fc_types["mat_undetected_pre_down"]
        if pre_fc >= hi_cutoff:
            return fc_types["mat_undetected_pre_up"]
        return fc_types["other"]
    if mat_fc <= low_cutoff and pre_fc >= hi_cutoff:
        return fc_types["mat_down_pre_up"]
    if mat_fc >= hi_cutoff and pre_fc >= higher_than_hi_factor * mat_fc:
        return fc_types["mat_up_pre_up"]
    if mat_fc <= low_cutoff:
        return fc_types["mat_down_only"]
    if mat_fc >= hi_cutoff:
        return fc_types["mat_up_only"]
    if pre_fc <= low_cutoff:
        return fc_types["pre_down_only"]
    if pre_fc >= hi_cutoff:
        return fc_types["pre_up_only"]
    return fc_types["other"]


def classify_fc_simple(mat_fc, pre_fc):
    change_type = classify_fc(mat_fc, pre_fc)
    if change_type in {fc_types["mat_down_pre_up"]}:
        return fc_types["mat_down_pre_up"]
    elif change_type in {
        fc_types["mat_down_pre_undetected"],
        fc_types["mat_down_only"],
    }:
        return "Only mature down"
    elif change_type in {
        fc_types["mat_undetected_pre_up"],
        fc_types["pre_up_only"],
        fc_types["mat_up_pre_up"],
    }:
        return "Only precursor up"
    else:
        return fc_types["other"]


data = pd.read_excel("extra/piRNA_summary.xlsx", sheet_name="piRNAs")
data.set_index("Name", inplace=True)
data["type"] = data[["Mat. FC (X2 RNAi/EV)", "Pre. FC (X2 RNAi/EV)"]].apply(
    lambda row: classify_fc_simple(
        row["Mat. FC (X2 RNAi/EV)"], row["Pre. FC (X2 RNAi/EV)"]
    ),
    axis="columns",
)
data["class"] = data.index.map(lambda name: "type2" if "type2" in name else "type1")
data_fill = data.fillna(0)

types_of_interest_strict = {
    # fc_types["mat_down_pre_undetected"],
    fc_types["mat_down_pre_up"],
    # fc_types["mat_up_pre_up"],
    # fc_types["mat_down_only"],
}

types_of_interest_loose = {
    fc_types["mat_down_pre_undetected"],
    fc_types["mat_down_pre_up"],
    # fc_types["mat_up_pre_up"],
    fc_types["mat_down_only"],
}

fig1 = plt.figure(dpi=600, figsize=(8, 8))
ax1 = fig1.gca()
sns.scatterplot(
    data=data_fill, x="Log2FC Mature", y="Log2FC Precursor", hue="type", ax=ax1, s=12
)
handles, labels = ax1.get_legend_handles_labels()
labels = [
    f"{label} = {data_fill[data_fill['type'] == label].count().min()}"
    for label in labels
]
handles, labels = zip(*sorted(zip(handles, labels), key=lambda x: x[1]))
ax1.legend(handles, labels, bbox_to_anchor=(1.025, 0.75))

log2_lc, log2_hc = np.log2(low_cutoff), np.log2(hi_cutoff)

ax1.axvline(x=0, linewidth=1, color="black")
ax1.axvline(x=log2_lc, linewidth=1, color="grey")
ax1.text(log2_lc - 0.55, ax1.get_ylim()[0] + 0.05, rf"{low_cutoff}$\times$")
ax1.axvline(x=log2_hc, linewidth=1, color="grey")
ax1.text(log2_hc + 0.1, ax1.get_ylim()[0] + 0.05, rf"{hi_cutoff}$\times$")
ax1.axhline(y=0, linewidth=1, color="black")
ax1.axhline(y=log2_lc, linewidth=1, color="grey")
ax1.text(ax1.get_xlim()[1] - 0.5, log2_lc - 0.25, rf"{low_cutoff}$\times$")
ax1.axhline(y=log2_hc, linewidth=1, color="grey")
ax1.text(ax1.get_xlim()[1] - 0.5, log2_hc + 0.15, rf"{hi_cutoff}$\times$")
ax1.set_xlabel(
    r"$\log_2\left(\mathrm{FC}\left(\frac{\mathrm{Mat.}_\mathrm{XRN2\:RNAi}}{\mathrm{Mat.}_\mathrm{EV}}\right)\right)$"
)
ax1.set_ylabel(
    r"$\log_2\left(\mathrm{FC}\left(\frac{\mathrm{Pre.}_\mathrm{XRN2\:RNAi}}{\mathrm{Pre.}_\mathrm{EV}}\right)\right)$"
)
fig1.savefig("extra/pirnas_mat_pre.svg", bbox_inches="tight")

fig2 = plt.figure(dpi=600, figsize=(15, 5))
ax2 = fig2.gca()
sub_data = data[data["Chr"] == "IV"].copy()
sub_data["x"] = sub_data["Start"] - sub_data["Start"].min()
sub_data.sort_values(by="x", inplace=True)

sns.scatterplot(data=sub_data, x="x", y="Log2FC Mature", hue="type", ax=ax2, s=12)

ymin, ymax = -5, 5
ax2.set_xlim(xmin=1e6)
ax2.set_ylim(ymin=ymin, ymax=ymax)
trans_y = lambda y: ((y - ymin) / (ymax - ymin))

for x, y, c in zip(sub_data["x"], sub_data["Log2FC Mature"], sub_data["type"]):
    if not pd.isna(y):
        if y < 0:
            ax2.axvline(
                x=x,
                ymin=trans_y(y),
                ymax=trans_y(0),
                color="blue" if c in types_of_interest_loose else "orange",
                linewidth=0.75 if c == fc_types["other"] else 1,
            )
        else:
            ax2.axvline(
                x=x,
                ymin=trans_y(0),
                ymax=trans_y(y),
                color="blue" if c in types_of_interest_loose else "orange",
                linewidth=0.75 if c == fc_types["other"] else 1,
            )
ax2.ticklabel_format(style="plain")
ax2.set_xlabel("Chr IV Coordinate")
ax2.set_ylabel(
    r"$\log_2\left(\mathrm{FC}\left(\frac{\mathrm{Mat.}_\mathrm{XRN2\:RNAi}}{\mathrm{Mat.}_\mathrm{EV}}\right)\right)$"
)
ax2.axhline(y=0, linewidth=1, color="black")
ax2.axhline(y=log2_hc, linewidth=1, color="grey")
ax2.axhline(y=log2_lc, linewidth=1, color="grey")
ax2.get_legend().remove()

fig2.savefig("extra/pirnas_chromosomal_distrib.svg", bbox_inches="tight")

all_detected_data = data[data["type"].isin(types_of_interest_strict)]
fig3 = plt.figure(dpi=600, figsize=(12, 3))
ax3, ax4 = fig3.subplots(nrows=1, ncols=2)
ax3 = sns.histplot(
    data=all_detected_data,
    x="EV Adj Len",
    binrange=(16, 28),
    discrete=True,
    stat="percent",
    ax=ax3,
)
# ax3.set_ylim(ymin=0.0, ymax=1.0)
ax3.set_title("EV")
ax3.set_xlabel("piRNA Species Length (nt)")
ax4 = sns.histplot(
    data=all_detected_data,
    x="X2 RNAi Adj Len",
    binrange=(16, 28),
    discrete=True,
    stat="percent",
    ax=ax4,
)
# ax4.set_ylim(ymin=0.0, ymax=1.0)
ax4.set_title("XRN2 RNAi")
ax4.set_xlabel("piRNA Species Length (nt)")

fig3.savefig("extra/pirnas_len_distrib.svg", bbox_inches="tight")

fig4 = plt.figure(figsize=(3,5), dpi=300)
ax4 = fig4.gca()
sns.boxplot(ax=ax4, data=all_detected_data[['EV Adj Len', 'X2 RNAi Adj Len']])
ax4.set_ylabel('Mean Read Length (nt)')
ax4.set_xticklabels(['EV', 'XRN2 RNAi'])
# ax4.set_title('Mean Read Length Distribution for XRN2 Target piRNAs')
fig4.savefig("extra/piRNAs_len_boxplot.svg", bbox_inches="tight")

def is_between(arr, low, high, low_incl=False, high_incl=True):
    if low_incl and high_incl:
        return (arr >= low) & (arr <= high)
    if low_incl:
        return (arr >= low) & (arr < high)
    if high_incl:
        return (arr > low) & (arr <= high)
    return (arr > low) & (arr < high)


type1_data = data[data["class"] == "type1"]
type2_data = data[data["class"] == "type2"]
print(
    "detected_mature_type1",
    type1_data[(type1_data["Mat. FC (X2 RNAi/EV)"] > 0)][
        "Mat. FC (X2 RNAi/EV)"
    ].count(),
)
print(
    "detected_mature_type2",
    type2_data[(type2_data["Mat. FC (X2 RNAi/EV)"] > 0)][
        "Mat. FC (X2 RNAi/EV)"
    ].count(),
)
print(
    "detected_pre_type1",
    type1_data[(type1_data["Pre. FC (X2 RNAi/EV)"] > 0)][
        "Pre. FC (X2 RNAi/EV)"
    ].count(),
)
print(
    "detected_pre_type2",
    type2_data[(type2_data["Pre. FC (X2 RNAi/EV)"] > 0)][
        "Pre. FC (X2 RNAi/EV)"
    ].count(),
)
print(
    "type1_mature_down",
    type1_data[is_between(type1_data["Mat. FC (X2 RNAi/EV)"], 0, low_cutoff)][
        "Mat. FC (X2 RNAi/EV)"
    ].count(),
)
print(
    "type2_mature_down",
    type2_data[is_between(type2_data["Mat. FC (X2 RNAi/EV)"], 0, low_cutoff)][
        "Mat. FC (X2 RNAi/EV)"
    ].count(),
)
print(
    "type1_pre_up",
    type1_data[type1_data["Pre. FC (X2 RNAi/EV)"] >= hi_cutoff][
        "Pre. FC (X2 RNAi/EV)"
    ].count(),
)
print(
    "type2_pre_up",
    type2_data[type2_data["Pre. FC (X2 RNAi/EV)"] >= hi_cutoff][
        "Pre. FC (X2 RNAi/EV)"
    ].count(),
)
print(
    "type1_mat_down_pre_up",
    type1_data[
        is_between(type1_data["Mat. FC (X2 RNAi/EV)"], 0, low_cutoff)
        & (type1_data["Pre. FC (X2 RNAi/EV)"] >= hi_cutoff)
    ]["Pre. FC (X2 RNAi/EV)"].count(),
)
print(
    "type2_mat_down_pre_up",
    type2_data[
        is_between(type2_data["Mat. FC (X2 RNAi/EV)"], 0, low_cutoff)
        & (type2_data["Pre. FC (X2 RNAi/EV)"] >= hi_cutoff)
    ]["Pre. FC (X2 RNAi/EV)"].count(),
)
print(
    data[data["type"] == fc_types["mat_down_pre_up"]]
    .sort_values(by="Pre. FC (X2 RNAi/EV)", ascending=False)
    .head(5)
)
# %%

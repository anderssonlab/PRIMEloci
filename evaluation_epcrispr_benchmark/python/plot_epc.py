import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def plot_count_zero(pos_df_sum_counts,
                    neg_df_sum_counts,
                    model, samples, colors,
                    P, N, output_pdf):
    # Masking values == 0
    pos_df_sum_counts_mask0 = pos_df_sum_counts[model][samples] == 0
    neg_df_sum_counts_mask0 = neg_df_sum_counts[model][samples] == 0

    # Counting zeros
    pos_count = pos_df_sum_counts_mask0.sum(axis=0)
    neg_count = neg_df_sum_counts_mask0.sum(axis=0)

    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    # Plot positive barplot
    sns.barplot(x=pos_count.index, y=pos_count.values, palette=colors, ax=ax1)
    ax1.set_title(f'Data that have COUNT = 0 \n EPC Positive : out of {P}')
    ax1.set_xlabel('Libraries')
    ax1.set_ylabel('Number')
    ax1.set_ylim(0, P)
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha="right")

    # Plot negative barplot
    sns.barplot(x=neg_count.index, y=neg_count.values, palette=colors, ax=ax2)
    ax2.set_title(f'Data that have COUNT = 0 \n EPC Negative : out of {N}')
    ax2.set_xlabel('Libraries')
    ax2.set_ylim(0, N)
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()

    # Save the figure as a PDF
    plt.savefig(output_pdf)

    # Show the plot
    plt.show()


def proba_confusion_matrix(df_pos, df_neg, threshold, P, N):

    predpos = df_pos >= threshold
    TP = predpos.sum(axis=0)
    FN = P - TP
    FN_cage = df_pos.isna().sum(axis=0)
    FN_pred = FN - FN_cage

    predneg = df_neg < threshold
    TN_pred = predneg.sum(axis=0)
    TN_cage = df_neg.isna().sum(axis=0)
    TN = TN_pred + TN_cage
    FP = N - TN

    confusion = pd.DataFrame({'Lib': TP.index,
                              'TP': TP,
                              'FN': FN,
                              'TN': TN,
                              'FP': FP,
                              'FN_cage': FN_cage,
                              'FN_pred': FN_pred,
                              'TN_cage': TN_cage,
                              'TN_pred': TN_pred})

    return confusion


def precaision_recall(confusion):
    precision = confusion['TP'] / (confusion['TP'] + confusion['FP'])
    recall = confusion['TP'] / (confusion['TP'] + confusion['FN'])
    df = pd.DataFrame({'precision': precision,
                       'recall': recall})

    return df


def plot_confusion(confusion, P, N, colors, save):

    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(16, 8))
    axes = axes.flatten()
    ylimit = [P, P, N, N, P, P, N, N]

    for i, column in enumerate(confusion.columns[1:]):
        # Exclude the first column ('Lib')
        sns.barplot(x='Lib', y=column,
                    data=confusion, ax=axes[i],
                    palette=colors)
        axes[i].set_title(f'{column}')
        axes[i].tick_params(axis='x', rotation=45)
        axes[i].set_ylim(0, ylimit[i])
    plt.tight_layout()
    if save is not False:
        plt.savefig(save + '.pdf')
    plt.show()


def wrap_cfs_pr_plot(df_pos, df_neg, threshole, P, N, colors, save=False):
    cfs = proba_confusion_matrix(df_pos, df_neg, threshole, P, N)
    pr = precaision_recall(cfs)
    plot_confusion(cfs, P, N, colors, save)
    return cfs, pr

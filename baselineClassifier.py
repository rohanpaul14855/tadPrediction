import os, csv, pickle
import matplotlib.pyplot as plt
from sklearn.metrics import f1_score, precision_score, recall_score, auc
from sklearn.metrics import accuracy_score, classification_report, roc_curve
from sklearn.metrics import precision_recall_curve, average_precision_score


# classifies a boundary as a TAD if it exists in at least one other cell line
# changing or to and
def main():
    with open('../data/pickle/pairwise_presence.unexpanded.pickle', 'rb') as f:
        df = pickle.load(f)
    for cell in df.columns.values:
        othercls = [k for k in df.columns.values if k != cell]
        # prediction = []
        # for i, row in df.iterrows():
        #     r = row[othercls]
        #     num_true = len([k for k in r if k == True])
        #     if num_true >= 3:
        #         prediction.append(True)
        #     else:
        #         prediction.append(False)
        prediction = df[othercls[0]] & df[othercls[1]] & df[othercls[2]] & df[othercls[3]]
        cr = classification_report(df[cell], prediction)
        # print(classification_report(df[cell], prediction),  '\n')
        # print(np.array([k.strip().split() for k in cr.split('\n') if len(k.strip()) > 0])); exit()
        # print('accuracy: {}\n'.format(accuracy_score(df[cell], prediction)))
        precision = precision_score(df[cell], prediction)
        recall = recall_score(df[cell], prediction)
        f1 = (2 * precision * recall) / (precision + recall)
        print('{}\t{}\t{}\t{}'.format(cell, precision, recall, f1))


def baseline():
    with open('../data/21/pickle/pairwise_presence.unexpanded.pickle', 'rb') as f:
        df = pickle.load(f)
    for cell in df.columns.values:
        if cell == 'BL':
            continue  # no data for BL

        # othercls = [k for k in df.columns.values if k != cell]
        # df2 = df[~(df[cell] & ~(df[othercls[0]] | df[othercls[1]] | df[othercls[2]] | df[othercls[3]]))] #True in cell and at least one other cell

        df2 = df.drop(cell, axis=1)
        any_true = [row.any() for _, row in df2.iterrows()]
        df2 = df2[any_true]  # drop all rows where every entry is False (not using cell line specific tads for ROC)
        dropped_idx = [df.ix[i].name for i, k in enumerate(any_true) if not k]
        cell_tad_presence = df[cell].drop(dropped_idx)  # real TADs of the cell line excluding cell line specific

        scores = [sum(row) for _, row in df2.iterrows()]  # get # of times a tad occurs in other cells
        fpr, tpr, threshold = roc_curve(cell_tad_presence, scores)
        precision, recall, _ = precision_recall_curve(cell_tad_presence, scores)
        roc_auc = auc(fpr, tpr)
        aupr = average_precision_score(cell_tad_presence, scores)
        plt.subplot(1, 2, 1)
        plt.plot(fpr, tpr, label='{} AUC: {:.2f}'.format(cell, round(roc_auc, 3)))
        plt.subplot(1, 2, 2)
        plt.plot(precision, recall, label='AUPR for {}: {}'.format(cell, round(aupr, 3)))

    # Plot ROC Curves
    plt.subplot(1, 2, 1)
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlabel('False Positive Rate')
    plt.xlabel('True Positive Rate')
    plt.title('Baseline AUROCs')

    # Plot PR curves
    plt.subplot(1, 2, 2)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='lower left')
    plt.title('Baseline AUPRs')

    plt.show()

    return


# main()
baseline()

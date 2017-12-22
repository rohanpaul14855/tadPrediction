import sklearn
from sklearn import linear_model
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle, os, csv
from collections import defaultdict
from sklearn.metrics import roc_auc_score, auc, roc_curve, average_precision_score, classification_report
from sklearn.metrics import precision_recall_curve, accuracy_score, f1_score, precision_score, recall_score
# from sklearn.model_selection import GridSearchCV

#Using binary features or not
BINARY = True


def readData(cell, binary=True):
    """
    :param cell: The cell line for which train and test data should be fetched
    :param binary: Whether to use binary or continuous features
    replace with style:[binary | continuous | average | aggregate]
    :return: returns the training and testing data as pandas dataframes
    """

    # path = '../data/main/unexpandedFeatures/train/{}'.format(cell)
    path = '../data/main/mainFeatures/{}'.format(cell)
    extension = '.bin.txt'
    train, test = [], []
    if not binary:
        extension = '.cont.txt'

    #training data
    for file in [k for k in os.listdir(path) if k.endswith(extension)]:
        fname = os.path.join(path, file)
        train.append(pd.read_csv(fname, index_col=0, header=0, sep='\t'))

    # Using bag of boundaries method for testing data
    # i.e. using same tads for training and testing
    # path = '../data/main/unexpandedFeatures/test/'
    path = '../data/main/BOBTestFeatures'
    for file in [k for k in os.listdir(path) if k.endswith(extension) and cell in k]:
        fname = os.path.join(path, file)
        test.append(pd.read_csv(fname, index_col=0, header=0, sep='\t'))

    train = pd.concat(train)
    test = pd.concat(test)

    ####
    # Add number of other cell lines each boundary is present in feature
    train_tads = train.index
    test_tads = test.index

    with open('../data/pickle/pairwise_presence.pickle', 'rb') as f:
        df = pickle.load(f)
        df = df.drop(cell, axis=1)

    train_counts = []
    for t in train_tads:
        numTrue = len([k for k in df.loc[t] if k == True])
        if numTrue == 0:
            exit(1)
        else:
            train_counts.append(numTrue)

    test_counts = []
    for t in test_tads:
        numTrue = len([k for k in df.loc[t] if k == True])
        if numTrue == 0:
            exit(1)
        else:
            test_counts.append(numTrue)
    train_counts = np.array(train_counts)
    test_counts = np.array(test_counts)
        # print(df.loc[t])
    train['counts'] = (train_counts - min(train_counts))/(max(train_counts) - min(train_counts))
    test['counts'] = (test_counts - min(test_counts))/(max(test_counts) - min(test_counts))

    #testing data
    # path = '../data/main/mainFeatures/testData/'
    # for file in [k for k in os.listdir(path) if k.endswith(extension) and cell in k]:
    #     fname = os.path.join(path, file)
    #     test.append(pd.read_csv(fname, index_col=0, header=0, sep='\t'))
    #
    # pospath = '../data/main/mainFeatures/testData/'
    # for file in [k for k in os.listdir(pospath) if k.endswith(extension) and cell in k and 'pos' in k]:
    #     fname = os.path.join(pospath, file)
    #     test.append(pd.read_csv(fname, index_col=0, header=0, sep='\t'))
    #
    # #One end in all tads is fixed with a real boundary
    # negpath = '../data/main/fixedNegativeFeatures'
    # for file in [k for k in os.listdir(negpath) if k.endswith(extension) and cell in k]:
    #     fname = os.path.join(negpath, file)
    #     test.append(pd.read_csv(fname, index_col=0, header=0, sep='\t'))

    return train, test




def plotCurves(fpr, tpr, cl, fn=None):
    if fn == None:
        fn = 'all features'
    roc_auc = auc(fpr, tpr)
    plt.title('ROC: Predicting {0} boundaries using {1}'.format(cl, fn))
    plt.plot(fpr, tpr, 'b', label='AUC = %0.2f' % roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')

    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    # plt.savefig('/Users/rohanpaul/Desktop/figs/pairwise_binary/{}_{}_cos.png'.format(cl, fn), bbox_inches='tight')
    plt.show()
    plt.close()


def classify(train, test, cl, binary=True):
    Xtrain = train[[k for k in train.columns.values if k!= 'Class']]
    ytrain = train['Class']
    Xtest = test[[k for k in test.columns.values if k!= 'Class']]
    ytest = test['Class']

    if not binary:  #normalize data
        Xtrain = (Xtrain - Xtrain.mean())/Xtrain.std()
        Xtest = (Xtest - Xtest.mean())/Xtest.std()

    bestAlpha = 0.01
    clf = linear_model.SGDClassifier(loss='modified_huber', class_weight='balanced', n_iter=10, n_jobs=-1,
                                     shuffle=True, alpha=bestAlpha)
    clf.fit(Xtrain, ytrain)
    scores = clf.predict_proba(Xtest)

    prediction = clf.predict(Xtest)


    posScores = scores[:, 1].reshape(-1, 1).tolist()
    weights = {}
    for i, feature in enumerate(Xtrain.columns.values):
        if feature == 'Correlation':
            weights[feature] = clf.coef_[0][i]
        elif feature == 'counts':
            weights[feature] = clf.coef_[0][i]
        else:
            weights[feature[:-2]] = clf.coef_[0][i]
    fpr, tpr, _ = roc_curve(ytest, posScores)
    precision, recall, _ = precision_recall_curve(ytest, posScores)
    print(cl)
    print('AUROC:', auc(fpr, tpr))
    print('AUPR:', average_precision_score(ytest, posScores))
    aupr = average_precision_score(ytest, posScores)
    print(classification_report(ytest, prediction))
    print()
    # print('Feature weights:  ')
    # for f, w in weights.items():
    #     print(f, w)
    # plotCurves(fpr, tpr, cl)
    # return weights
    return weights, fpr, tpr, precision, recall, aupr



def SVM(train, test, cl, binary=True):

    Xtrain = train[[k for k in train.columns.values if k != 'Class']]
    ytrain = train['Class']
    Xtest = test[[k for k in test.columns.values if k != 'Class']]
    ytest = test['Class']

    if not binary:  #normalize data
        Xtrain = (Xtrain - Xtrain.mean())/Xtrain.std()
        Xtest = (Xtest - Xtest.mean())/Xtest.std()

    clf = sklearn.svm.LinearSVC(class_weight='balanced', C=1)
    clf.fit(Xtrain, ytrain)
    prediction = clf.predict(Xtest)
    print(cl)
    print(classification_report(ytest, prediction))
    print('\n')


def main():
    weights = defaultdict(list)
    for cell in ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']:
        train, test = readData(cell, binary=BINARY)
        w, fpr, tpr, precision, recall, aupr = classify(train, test, cell, binary=BINARY)
        for t in w:
            weights[t].append(w[t])
        plt.subplot(1, 2, 1)
        plt.plot(fpr, tpr, label='AUROC for {}: {}'.format(cell, round(auc(fpr, tpr), 3)))
        plt.subplot(1, 2, 2)
        plt.plot(precision, recall, label='AUPR for {}: {}'.format(cell, round(aupr, 3)))
    # plot ROC curves
    plt.subplot(1,2,1)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')
    plt.title('Bag of Boundaries ROCs')
    plt.plot([0, 1], [0, 1], 'k--')
    # plot PR curves
    plt.subplot(1, 2, 2)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='lower left')
    plt.title('Bag of Boundaries AUPRs')
    # plt.plot([0, 1], [1, 0], 'k--')
    plt.show()

    # plot feature weights
    names = list(weights.keys())
    ind = np.arange(len(names))
    width = 0.35
    p1 = plt.bar(ind, [np.mean(a) for a in weights.values()], width, color='0.5')
    plt.ylabel('Feature Weights')
    plt.title('Feature Importance')
    plt.xticks(ind, names, rotation='vertical')
    plt.show()
    plt.close()

main()

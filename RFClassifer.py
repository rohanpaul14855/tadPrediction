import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np
from sklearn.metrics import auc, precision_recall_curve, average_precision_score
from sklearn.metrics import roc_curve, f1_score, classification_report
from sklearn.metrics import precision_score, recall_score
from collections import defaultdict
import seaborn as sns
import csv, os
import pickle
from sklearn import cross_validation

BINARY = True

def plotCurves(fpr, tpr, cl, fn=None):
    if len(fpr) != 0:
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, 'b', label='AUC = %0.2f' % roc_auc)
    if fn == None:
        fn = 'all features'
    plt.title('ROC: Predicting {0} boundaries using {1}'.format(cl, fn))
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    ax = plt.gca()
    ax.grid(True)
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    # plt.savefig('/Users/rohanpaul/Desktop/figs/pairwise_binary/{}_{}_cos.png'.format(cl, fn), bbox_inches='tight')
    plt.show()
    plt.close()


def OneVsAll():
    w = defaultdict(list)
    cells = ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']
    for cell in cells:
        toconcat = []

        ocl = [k for k in cells if k != cell]
        for t in ocl:
            path = 'PairwiseFeatures(global)/{}_testData'.format(t)
            df = pd.read_csv(path, header=0, index_col=0, sep='\t')
            toconcat.append(df)
        train = pd.concat(toconcat)
        test = pd.read_csv('PairwiseFeatures(global)/{}_testData'.format(cell), header=0,
                           index_col=0, sep='\t')

        # parameters = {'n_estimators': [10, 100, 1000],
        #               'criterion': ['gini', 'entropy']}
        # rf = RandomForestClassifier(n_jobs=-1)
        # clf = GridSearchCV(rf, parameters)

        clf = RandomForestClassifier(n_estimators=200, criterion='gini',
                                     n_jobs=-1)

        Xtrain = train[train.columns.values[:-1]]
        ytrain = train['Class']
        Xtest = test[test.columns.values[:-1]]
        ytest = test['Class']

        clf.fit(Xtrain, ytrain)
        # print(clf.best_params_)
        scores = clf.predict_proba(Xtest)
        posScores = scores[:, 1].reshape(scores.shape[0])
        fpr, tpr, _ = roc_curve(ytest, posScores)


        ##############PLOT################

        roc_auc = round(auc(fpr, tpr), 3)
        print('AUROC for {}: {}'.format(cell, roc_auc))
        plt.plot(fpr, tpr, label='AUC for {} = {}'.format(roc_auc, cell))
    plt.title('One vs. All'.format(cell))
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    # plt.savefig('/Users/rohanpaul/Desktop/OnevAll.png'.format(cell), bbox_inches='tight')
    plt.show()
    plt.close()


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


def singleCellLine():
    w = defaultdict(list)
    for cl in ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']:
        master = []
        with open('../PairwiseFeatures(global)/{}_testData'.format(cl), 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for line in reader:
                    master.append(list(line))
        header = master.pop(0)
        # master = master[1:]
        master = pd.DataFrame(master, columns=header)
        # master.set_index('featureID')
        y = master['Class']
        X = master.drop('Class', axis=1)
        cv = cross_validation.StratifiedKFold(y, n_folds=5, shuffle=True)
        y = y.values; X = X.drop('featureID', axis=1).values
        auroc = 0
        # print(X[:10]); exit()
        plt.plot([0, 1], [0, 1], 'r--')
        plt.xlabel('False Positive Rate')
        plt.title('{} Only ROC Curve'.format(cl))
        plt.ylabel('True Positive Rate')

        for fold, (train_index, test_index) in enumerate(cv):
            y_train, y_test = y[train_index], y[test_index]
            x_train, x_test = X[train_index], X[test_index]

            clf = RandomForestClassifier(n_estimators=200, criterion='gini',
                                         n_jobs=-1, max_depth=5)

            clf.fit(x_train, y_train)
            scores = clf.predict_proba(x_test)
            posScores = scores[:, 1].reshape(-1, 1).tolist()
            # for i, name in enumerate(master.columns.values[1:-1]):
            #     if name != 'Correlation':
            #         w[name[:-2]].append(clf.coef_[0][i])
            #     else:
            #         w[name].append(clf.coef_[0][i])
            fpr, tpr, _ = roc_curve(y_test, posScores, pos_label='1')
            auc_roc = round(auc(fpr, tpr), 3)
            plt.plot(fpr, tpr, label='ROC on fold {}: {}'.format(fold, auc_roc))
            print('AUC on fold {} for {}:'.format(fold, cl), auc_roc)
            auroc += auc(fpr, tpr)
        # plotCurves([], [], cl)
        plt.legend(loc='lower right')

        plt.show()
        print('Average AUC for {} = {}'.format(cl, auroc/5))

    # avw = []
    # for key, value in w.items():
    #     avw.append(sum(value)/len(value))
    # names = list(w.keys())
    # ind = np.arange(len(avw))
    # width = 0.35
    # p1 = plt.bar(ind, avw, width, color='0.5')
    # plt.ylabel('Feature Weights')
    # plt.title('Feature Importance')
    # plt.xticks(ind, names, rotation='vertical')
    # plt.show()
    # plt.close()


def classify(train, test, cl, binary):

    Xtrain = train[[k for k in train.columns.values if k != 'Class']]
    ytrain = train['Class']
    Xtest = test[[k for k in test.columns.values if k != 'Class']]
    ytest = test['Class']

    if not binary:  #normalize data
        Xtrain = (Xtrain - Xtrain.mean())/(Xtrain.max() - Xtrain.min())
        Xtest = (Xtest - Xtest.mean())/(Xtest.max() - Xtest.min())


    clf = RandomForestClassifier(n_estimators=200, criterion='gini',
                                 n_jobs=-1, class_weight='balanced',
                                 max_depth=6)
    clf.fit(Xtrain, ytrain)
    scores = clf.predict_proba(Xtest)
    prediction = clf.predict(Xtest)
    precision = precision_score(ytest, prediction)
    recall = recall_score(ytest, prediction)
    f1 = f1_score(ytest, prediction)

    print('{}\t{}\t{}\t{}'.format(cl, precision, recall, f1))

    posScores = scores[:, 1].reshape(-1, 1).tolist()
    weights = {}
    # for i, feature in enumerate(Xtrain.columns.values):
    #     if feature == 'Correlation':
    #         weights[feature] = clf.coef_[0][i]
    #     else:
    #         weights[feature[:-2]] = clf.coef_[0][i]
    fpr, tpr, _ = roc_curve(ytest, posScores)
    precision, recall, _ = precision_recall_curve(ytest, posScores)
    # print(cl)
    print('AUROC:', auc(fpr, tpr))
    aupr = average_precision_score(ytest, posScores)
    # print(classification_report(ytest, prediction))
    # print()
    # print('Feature weights:  ')
    # for f, w in weights.items():
    #     print(f, w)
    # plotCurves(fpr, tpr, cl)
    return fpr, tpr, precision, recall, aupr


def read21Data(cell):
    # path = '../data/main/unexpandedFeatures/train/{}'.format(cell)
    path = '../data/21/features/train/{}'.format(cell)
    train, test = [], []

    for subdir in os.listdir(path):
        for feature_file in os.listdir(os.path.join(path, subdir)):
            if feature_file.endswith('.features'):
                fname = os.path.join(path, subdir, feature_file)
                train.append(pd.read_csv(fname, index_col=0, header=0, sep='\t'))

    path = '../data/21/features/test/'
    for feature_file in os.listdir(path):
        if cell in feature_file and feature_file.endswith('.features'):
            fname = os.path.join(path, feature_file)
            test.append(pd.read_csv(fname, index_col=0, header=0, sep='\t'))

    train = pd.concat(train)
    test = pd.concat(test)

    ####
    # Add number of other cell lines each boundary is present in feature
    train_tads = train.index
    test_tads = test.index

    with open('../data/21/pickle/pairwise_presence.unexpandedcont.pickle', 'rb') as f:
        df = pickle.load(f)
        df = df.drop(cell, axis=1)

    train_counts = []
    for t in train_tads:
        numTrue = len([k for k in df.loc[t] if k == True])
        if numTrue == 0:
            print('here1')
            exit(1)
        else:
            train_counts.append(numTrue)

    test_counts = []
    for t in test_tads:
        numTrue = len([k for k in df.loc[t] if k == True])
        if numTrue == 0:
            print('here2')
            exit(1)
        else:
            test_counts.append(numTrue)

    train_counts = np.array(train_counts)
    test_counts = np.array(test_counts)
    train['counts'] = (train_counts - min(train_counts)) / (max(train_counts) - min(train_counts))
    test['counts'] = (test_counts - min(test_counts)) / (max(test_counts) - min(test_counts))

    # testing data
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


def main():
    weights = defaultdict(list)
    cells = 'AD AO CO GM12878 H1 HC IMR90 LG LI LV MES MSC NPC OV PA PO RV SB SX TRO'.split()
    for cell in cells:
        train, test = read21Data(cell)
        fpr, tpr, precision, recall, aupr = classify(train, test, cell, binary=BINARY)
        plt.subplot(1, 2, 1)
        plt.plot(fpr, tpr, label='ROC for {}: {}'.format(cell, round(auc(fpr, tpr), 2)))
        plt.subplot(1, 2, 2)
        plt.plot(precision, recall, label='AUPR for {}: {}'.format(cell, round(aupr, 3)))
        # plot ROC curves
    plt.subplot(1, 2, 1)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')
    plt.title('Bag of Boundaries ROCs')
    plt.plot([0, 1], [0, 1], 'k--')
    # plot PR curves
    plt.subplot(1, 2, 2)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='upper left')
    plt.title('Bag of Boundaries AUPRs')
    plt.show()


main()
# singleCellLine()
# OneVsAll()



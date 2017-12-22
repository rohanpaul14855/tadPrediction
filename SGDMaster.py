#############################################################################
#@author = Rohan Paul
#07/26/2016
#Plots ROC curves for each feature
#############################################################################
from collections import defaultdict
import numpy as np
from sklearn import preprocessing
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, auc, roc_curve, average_precision_score
from sklearn.grid_search import GridSearchCV
from sklearn.cross_validation import train_test_split
from sklearn.model_selection import train_test_split
from sklearn import linear_model
from sklearn import cross_validation
import csv, sys, os
import pickle
import subprocess


######Calls cryptic pairs and stores 5C style pairs in a file
######Returns file paths
def callCryptic(pos_index, neg_index, cl):
    basePath = '/Users/rohanpaul/Dropbox/AyLab/'
    for i, name in zip([pos_index, neg_index], ['pos', 'neg']):
        Pairs5C = []
        for line in i:
            line = line.split('-')
            bound1 = line[0].split('_')
            bound2 = line[1].split('_')
            cb1 = '5C_000_ENm000_FOR_000|hg19|{}:{}-{}'.format(bound1[0], bound1[1], bound1[2])
            cb2 = '5C_000_ENm000_REV_000|hg19|{}:{}-{}'.format(bound2[0], bound2[1], bound2[2])
            Pairs5C.append([cb1, cb2])
        outfile = basePath + 'LocalPairsCryptic/{}_{}.txt'.format(cl, name)
        writer = csv.writer(open(outfile, 'w'), delimiter='\t')
        writer.writerows(Pairs5C)
    return (basePath + 'LocalPairsCryptic/{}_{}.txt'.format(cl, 'pos'),
            basePath + 'LocalPairsCryptic/{}_{}.txt'.format(cl, 'neg'))


###### Takes a list pos, a list neg
###### Calls Cryptic Pairs.py
###### Calls C++ program to generate features for
###### Returns path of feature files
def genFeatures(pos, neg, cl):
    paths = []
    os.chdir('/Users/rohanpaul/Dropbox/AyLab/forwebsite/code/generatefeatures')
    basePath = '/Users/rohanpaul/Dropbox/AyLab/'
    bashCommand = ['./genFeatures', 'pairFile', 'dataLocationFile', 'output', 'label', 'no', 'concat',
                   'null', 'yes', 'binary']
    #set pairFile location
    pairFile = callCryptic(pos.index, neg.index, cl)
    for name, varname, label in zip([pos, neg], ['pos', 'neg'], ['1', '0']):
        bashCommand[1] = pairFile[0]
        #set reg element file location
        bashCommand[2] = basePath + 'locs/{}_locs'.format(cl)
        bashCommand[3] = basePath + 'localOut/{}_positive_expanded.txt'.format(cl)
        bashCommand[4] = label
        execCommand = open('tempfiletoexec.bash', 'w')
        execCommand.write('#!/usr/bin/env bash')
        execCommand.write('\n')
        execCommand.write(' '.join(bashCommand))
        if not os.path.isfile(basePath + 'localPairFeatures/{}_{}_localPairFeatures.txt'.format(cl, varname)):
            ret = subprocess.check_output(['bash', 'tempfiletoexec.bash'])#, stdout=subprocess.PIPE)
            print('ret', '\n', ret)
        else:
            print('Already run before, check localPairFeatures/ directory and rerun if needed')
        execCommand.close()
    paths.append(basePath + 'localPairFeatures/{}_{}_localPairFeatures.txt'.format(cl, 'pos'))
    paths.append(basePath + 'localPairFeatures/{}_{}_localPairFeatures.txt'.format(cl, 'neg'))
    os.chdir(basePath)
    return(paths)


###############################################################################
#   Splits data by cell line for training and testing                         #
#   returns Xtrain, Xtest, ytrain, ytest split by cellline                    #
###############################################################################
def returnTrainTest_local(cl):
    cells = ['HeLa', 'HMEC', 'HUVEC', 'K562', 'NHEK']
    all_features = ['ctcf_E', 'h2afz_E', 'h3k27ac_E', 'h3k27me3_E', 'h3k36me3_E', 'h3k4me1_E',
                    'h3k4me2_E', 'h3k4me3_E', 'h3k79me2_E', 'h3k9ac_E', 'h3k9me3_E', 'h4k20me1_E', 'ctcf_P',
                    'h2afz_P', 'h3k27ac_P', 'h3k27me3_P', 'h3k36me3_P', 'h3k4me1_P', 'h3k4me2_P', 'h3k4me3_P',
                    'h3k79me2_P', 'h3k9ac_P', 'h3k9me3_P', 'h4k20me1_P', 'Correlation']

    with open('pickle/pairwise_presence.pickle', 'rb') as f:
        master = pickle.load(f)

    pos = master[master[cl] == True]
    pos_boundaries = list(map(lambda x: cl + '_' + x, pos.index))
    neg = master[master[cl] == False]

    for i in range(neg.shape[0] - 1):
        trues = neg.columns[neg.iloc[i] == True]


    # clIdx = [i+1 for i, k in enumerate(cells) if k == cl][0]
    # for line in master:
    #     numTrue = len([k for k in line if k == True])
    #     if line[clIdx] == True:
    #         pos.append('{}_{}'.format(cl, line[0]))
    #     elif line[clIdx] == False and numTrue >= 1:
    #         idx = [j for j, k in enumerate(line) if k == True][0]
    #         neg.append('{}_{}'.format(cells[idx-1], line[0]))
    featurePaths = genFeatures(pos, neg, cl)

    posMaster = pd.read_csv(featurePaths[0], sep='\t')
    posMaster = posMaster.convert_objects(convert_numeric=True, copy=True)
    posLabels = np.ones(posMaster.shape[0], dtype=int)

    negMaster = pd.read_csv(featurePaths[1], sep='\t')
    negMaster = negMaster.convert_objects(convert_numeric=True, copy=True)
    negLabels = np.zeroes(negMaster.shape[0], dtype=int)

    pos_x_train, pos_x_test, pos_y_train, pos_y_test = train_test_split(posMaster, posLabels, test_size=0.2)
    neg_x_train, neg_x_test, neg_y_train, neg_y_test = train_test_split(negMaster, negLabels, test_size=0.2)

    X_train = pd.concat([pos_x_train, neg_x_train])
    y_train = pd.concat([pd.Series(pos_y_train), pd.Series(neg_y_train)])
    X_test = pd.concat([pos_x_test, neg_x_test])
    y_test = pd.concat([pd.Series(pos_y_test), pd.Series(neg_y_test)])

    return X_train, X_test, y_train, y_test


def returnTrainTest_pairs(cl, feature=None, correlation=None):
    all_features = ['featureID', 'ctcf_E', 'h2afz_E', 'h3k27ac_E', 'h3k27me3_E', 'h3k36me3_E', 'h3k4me1_E',
                    'h3k4me2_E', 'h3k4me3_E', 'h3k79me2_E', 'h3k9ac_E', 'h3k9me3_E', 'h4k20me1_E', 'ctcf_P',
                    'h2afz_P', 'h3k27ac_P', 'h3k27me3_P', 'h3k36me3_P', 'h3k4me1_P', 'h3k4me2_P', 'h3k4me3_P',
                    'h3k79me2_P', 'h3k9ac_P', 'h3k9me3_P', 'h4k20me1_P', 'Correlation', 'Cosine']
    master = pd.read_csv('PairwiseFeatures(global)/master_pairwise_cos.txt', sep='\t')
    master = master.convert_objects(convert_numeric=True, copy=True)

    train = master[~master['featureID'].str.contains(cl)]
    test = master[master['featureID'].str.contains(cl)]

    Xtrain = train[all_features]
    Xtrain.set_index('featureID', inplace=True)
    ytrain = train['Class']
    Xtest = test[all_features]
    Xtest.set_index('featureID', inplace=True)
    ytest = test['Class']

    if correlation == None or correlation == 'Pearson':
        Xtrain.drop('Cosine', axis=1, inplace=True)
        Xtest.drop('Cosine', axis=1, inplace=True)
    elif correlation == 'Cosine':
        Xtrain.drop('Correlation', axis=1, inplace=True)
        Xtest.drop('Correlation', axis=1, inplace=True)
    else:
        sys.exit('Enter either \'Pearson\' or \'Cosine\' for correlation argument')
    return Xtrain, Xtest, ytrain, ytest

    # test = pd.read_csv('PairwiseFeatures(global)/{}_testData_cos'.format(cl), sep='\t')
    # test = test.convert_objects(convert_numeric=True, copy=True)
    # trainMaster, testMaster = [], []
    # with open('PairwiseFeatures(global)/master_pairwise', 'r') as trainData:
    # train = csv.reader(trainData, delimiter='\t')
    # trainMaster = list(train)
    # for i in pd.isnull(train).values:
    #     for j in i:
    #         # if j == 'False':
    #         #     print('garbl')
    #         if j == True:
    #             print(i,)
    # with open('PairwiseFeatures(global)/HeLa_testData', 'r') as testData:
    # test = csv.reader(testData, delimiter='\t')
    # testMaster = list(test)

    # for line in trainMaster:
    #     for i, element in enumerate(line):
    #         try:
    #             line[i] = float(element)
    #         except:
    #             continue
    # for line in testMaster:
    #     for i, element in enumerate(line):
    #         try:
    #             line[i] = float(element)
    #         except:
    #             continue

    # if feature != None:
    #     pass
        #implement this part for making classifiers for a single feature
    # else:
        # header = trainMaster.pop(0)
        # trainMaster = pd.DataFrame(trainMaster, columns=header)
        # testMaster = pd.DataFrame(testMaster, columns=header)

    # if correlation == None or correlation == 'Pearson':
    #     train.drop('Cosine', inplace=True, axis=1)
    #     test.drop('Cosine', inplace=True, axis=1)
    # elif correlation == 'Cosine':
    #     train.drop('Correlation', axis=1, inplace=True)
    #     test.drop('Correlation', axis=1, inplace=True)
    # else:
    #     sys.exit('Enter either \'Pearson\' or \'Cosine\' for correlation argument')


def returnTrainTest_single(cl, feature = None):
    trainMaster = []; testMaster = []

    with open('featureSelection/master_normed_DNase.txt','r') as f:
        master = csv.reader(f, delimiter='\t')
        for row in master:
            if cl in row[0]:
                testMaster.append(row)
            else:
                trainMaster.append(row)
    for row in trainMaster:
        for i, elem in enumerate(row):
            try:
                row[i] = float(elem)
            except:
                continue

    for row in testMaster:
        for i, elem in enumerate(row):
            try:
                row[i] = float(elem)
            except:
                continue

    if feature != None:
        fIdx = trainMaster[0].index(feature)
        Xtrain = [row[fIdx] for row in trainMaster][1:]
        Xtrain = np.asarray(Xtrain)
        Xtrain = Xtrain.reshape(-1,1)       #reshaping for singular feature

        ytrain = [row[1] for row in trainMaster][1:]
        ytrain = np.asarray(ytrain)
        ytrain.reshape(-1,1)

        Xtest = [row[fIdx] for row in testMaster][1:]
        Xtest = np.asarray(Xtest)
        Xtest = Xtest.reshape(-1, 1)

        ytest = [row[1] for row in testMaster][1:]
        ytest = np.asarray(ytest)
        ytest = ytest.reshape(-1,1)

    else:
        header = trainMaster.pop(0)

        trainMaster = pd.DataFrame(trainMaster, columns=header)
        # idx = trainMaster.drop('featureID', axis=1)
        #trainMaster.drop('featureID',1, inplace=True)
        #trainMaster.set_index('featureID', inplace=True)
        #######################################################################################################################
        # X_train = trainMaster


        testMaster = pd.DataFrame(testMaster, columns=header)
        # idx = testMaster.drop('featureID',axis=1)
        #testMaster.drop('featureID', 1, inplace=True)
        #testMaster.set_index('featureID', inplace=True)

        highROCs= ['featureID', 'ctcfAggSignal', 'ctcfPeakCounts',
                   'ctcfSignalStrength', 'h2afzAggSignal', 'h2afzPeakCounts',
                   'h2afzSignalStrength', 'h3k27acAggSignal' ,'h3k27acPeakCounts',
                   'h3k27acSignalStrength' , 'h3k36me3PeakCounts', 'h3k36me3SignalStrength',
                   'h3k36meAggSignal', 'h3k4me1AggSignal','h3k4me1PeakCounts',
                   'h3k4me1SignalStrength', 'h3k4me2AggSignal', 'h3k4me2PeakCounts',
                   'h3k4me2SignalStrength', 'h3k4me3AggSignal', 'h3k4me3PeakCounts',
                   'h3k4me3SignalStrength', 'h3k79me2AggSignal', 'h3k79me2PeakCounts',
                   'h3k79me2SignalStrength', 'h3k9acAggSignal', 'h3k9acPeakCounts',
                   'h3k9acSignalStrength', 'h3k9me3AggSignal',
                     'h4k20me1AggSignal' ,
                   'h4k20me1PeakCounts','h4k20me1SignalStrength', 'DNaseSignal']
        #pft = ['ctcfAggSignal', 'ctcfSignalStrength', 'h3k9acAggSignal', 'h2afzAggSignal', 'ctcfPeakCounts',
        #           'h3k27acAggSignal', 'h3k4me1AggSignal', 'h3k4me2AggSignal', 'h3k4me3AggSignal',
        #                            'h3k79me2AggSignal', 'h4k20me1AggSignal']


        Xtrain = trainMaster[highROCs]
        Xtrain.set_index('featureID', inplace=True)
        ytrain = trainMaster['Class']
        Xtest = testMaster[highROCs]
        Xtest.set_index('featureID', inplace=True)
        ytest = testMaster['Class']

    return Xtrain, Xtest, ytrain, ytest


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



def readFeatures(cl):
    posMaster = pd.read_csv('{}_positiveFeatures.txt'.format(cl), sep='\t', header=0, index_col=0)
    posLabels = posMaster['Class']
    posMaster.drop('Class', axis=1, inplace=True)
    posMaster = posMaster.convert_objects(convert_numeric=True, copy=True)


    negMaster = pd.read_csv('{}_negativeFeatures.txt'.format(cl), sep='\t', header=0, index_col=0)
    negLabels = negMaster['Class']
    negMaster.drop('Class', axis=1, inplace=True)
    negMaster = negMaster.convert_objects(convert_numeric=True, copy=True)

    pos_x_train, pos_x_test, pos_y_train, pos_y_test = train_test_split(posMaster, posLabels, test_size=0.2)
    neg_x_train, neg_x_test, neg_y_train, neg_y_test = train_test_split(negMaster, negLabels, test_size=0.2)

    X_train = pd.concat([pos_x_train, neg_x_train])
    y_train = pd.concat([pd.Series(pos_y_train), pd.Series(neg_y_train)])
    X_test = pd.concat([pos_x_test, neg_x_test])
    y_test = pd.concat([pd.Series(pos_y_test), pd.Series(neg_y_test)])

    return X_train, X_test, y_train, y_test

def main():

    for cl in ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']:
    # for cl in ['HeLa']:
        # Xtrain, Xtest, ytrain, ytest = readFeatures(cl)
        posTrain = pd.read_csv('trainData/{}_pos_features.txt'.format(cl), sep='\t', index_col=0, header=0)
        negTrain = pd.read_csv('trainData/{}_neg_features.txt'.format(cl), sep='\t', index_col=0, header=0)


        posLabels = posTrain['Class']
        posTrain.drop('Class', axis=1, inplace=True)
        negLabels = negTrain['Class']
        negTrain.drop('Class', axis=1, inplace=True)

        pos_x_train, pos_x_test, pos_y_train, pos_y_test = train_test_split(posTrain, posLabels, test_size=0.2)
        neg_x_train, neg_x_test, neg_y_train, neg_y_test = train_test_split(negTrain, negLabels, test_size=0.2)

        Xtrain = pd.concat([pos_x_train, neg_x_train])
        ytrain = pd.concat([pd.Series(pos_y_train), pd.Series(neg_y_train)])
        Xtest = pd.concat([pos_x_test, neg_x_test])
        ytest = pd.concat([pd.Series(pos_y_test), pd.Series(neg_y_test)])

        bestAlpha = 0.01
        clf = linear_model.SGDClassifier(loss='modified_huber', class_weight='balanced', n_iter=10, n_jobs=-1,
                                         shuffle=True, alpha=bestAlpha)

        clf.fit(Xtrain, ytrain)
        scores = clf.predict_proba(Xtest)
        # pred = clf.predict(Xtrain)
        # print(sum(pred != ytrain)); exit()
        # good = []
        # bad = []
        #
        # for i, prediction in enumerate(pred):
        #     if prediction != ytest[i]:
        #         bad.append(Xtest.iloc[i])
        #     else:
        #         good.append(Xtest.iloc[i])
        #
        # goodPredictions = pd.DataFrame(good)
        # badPredictions = pd.DataFrame(bad)
        #
        # goodPredictions.to_csv('goodPreds.csv', sep='\t')
        # badPredictions.to_csv('badPreds.csv', sep='\t')
        # TP, FP, TN, FN = [], [], [], []
        # for i, prediction in enumerate(pred):
        #     if prediction == -1 and ytest[i] == -1:
        #         TN.append(Xtest.iloc[i])
        #     elif prediction == 1 and ytest[i] == 1:
        #         TP.append(Xtest.iloc[i])
        #     elif prediction == -1 and ytest[i] == 1:
        #         FN.append(Xtest.iloc[i])
        #     elif prediction == 1 and ytest[i] == -1:
        #         FP.append(Xtest.iloc[i])
        # pd.DataFrame(TP).to_csv('TP.csv', sep='\t')
        # pd.DataFrame(FP).to_csv('FP.csv', sep='\t')
        # pd.DataFrame(TN).to_csv('TN.csv', sep='\t')
        # pd.DataFrame(FN).to_csv('FN.csv', sep='\t')
        #
        posScores = scores[:, 1].reshape(-1, 1).tolist()

        fpr, tpr, _ = roc_curve(ytest, posScores)
        print('auc:', auc(fpr, tpr))
        print('pr:', average_precision_score(ytest, posScores))
        plotCurves(fpr, tpr, cl)


def singleCellLine(binary=True):
    w = defaultdict(list)
    extension = 'bin.txt'
    if not binary:
        extension = 'cont.txt'

    files = [k for k in os.listdir('../data/single/singleFeatures/') if k.endswith(extension)]
    for cl in ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']:
        cl_files = [os.path.join('../data/single/singleFeatures', k) for k in files if cl in k]
        data = []
        for file in cl_files:
            df = pd.read_csv(file, sep='\t', index_col=0, header=0)
            data.append(df)
        master = pd.concat(data)


        y = master['Class']
        X = master.drop('Class', axis=1)
        X = (X - X.mean())/X.std()
        cv = cross_validation.StratifiedKFold(y, n_folds=5, shuffle=True)
        auroc = 0
        '''plt.title('Predicting {0} boundaries using all features'.format(cl))
        plt.legend(loc='lower right')
        plt.plot([0, 1], [0, 1], 'r--')
        plt.ylabel('True Positive Rate')
        plt.xlabel('False Positive Rate')'''

        for fold, (train_index, test_index) in enumerate(cv):

            y_train, y_test = y.iloc[train_index], y.iloc[test_index]
            x_train, x_test = X.iloc[train_index], X.iloc[test_index]
            clf = linear_model.SGDClassifier(loss='log', class_weight='balanced', n_iter=10, n_jobs=-1,
                                             shuffle=True, alpha=1e-2)
            clf.fit(x_train, y_train)
            scores = clf.predict_proba(x_test)
            posScores = scores[:, 1].reshape(-1, 1).tolist()

            for i, name in enumerate(master.columns.values[1:-1]):
                if name != 'Correlation':
                    w[name[:-2]].append(clf.coef_[0][i])
                else:
                    w[name].append(clf.coef_[0][i])
            fpr, tpr, _ = roc_curve(y_test, posScores, pos_label=1)
            print('AUC on fold {} for {}:'.format(fold, cl), auc(fpr, tpr))
            auroc += auc(fpr, tpr)
            # plotCurves(fpr, tpr, cl)
        print('Average AUC for {} = {}'.format(cl, auroc/5))
    avw = []
    for key, value in w.items():
        avw.append(sum(value)/len(value))
    names = list(w.keys())
    ind = np.arange(len(avw))
    width = 0.35
    p1 = plt.bar(ind, avw, width, color='0.5')
    plt.ylabel('Feature Weights')
    plt.title('Feature Importance')
    plt.xticks(ind, names, rotation='vertical')
    plt.show()
    plt.close()
    # plt.legend()

def OneVsAll():
    w = defaultdict(list)
    cells = ['HeLa', 'HMEC', 'HUVEC', 'NHEK', 'K562']
    files = [k for k in os.listdir('../data/single/singleFeatures') if k.endswith('bin.txt')]
    for cell in cells:
        toconcat = []
        ocl = [k for k in files if cell not in k]
        for t in ocl:
            path = os.path.join('../data/single/singleFeatures/', t)
            df = pd.read_csv(path, header=0, index_col=0, sep='\t')
            toconcat.append(df)
        train = pd.concat(toconcat)
        test_cls = [k for k in files if cell in k]
        toconcat = []
        for t in test_cls:
            path = os.path.join('../data/single/singleFeatures/', t)
            df = pd.read_csv(path, header=0, index_col=0, sep='\t')
            toconcat.append(df)
        test = pd.concat(toconcat)

        clf = linear_model.SGDClassifier(loss='modified_huber', class_weight='balanced', n_iter=10,
                                         n_jobs=-1, shuffle=True, alpha=1e-2)
        x_train = train[train.columns.values[:-1]]
        y_train = train['Class']
        x_test = test[test.columns.values[:-1]]
        y_test = test['Class']
        clf.fit(x_train, y_train)
        scores = clf.predict_proba(x_test)
        prediction = clf.predict(x_test)
        posScores = scores[:, 1].reshape(scores.shape[0])
        fpr, tpr, _ = roc_curve(y_test, posScores)
        # plotCurves(fpr, tpr, cell)
        ##############PLOT################
        roc_auc = auc(fpr, tpr)
        print('AUROC for {}: {}'.format(cell, roc_auc))
        plt.plot(fpr, tpr, label='AUC = %0.2f' % roc_auc)
        plt.title('One vs. All: {}'.format(cell))
        plt.legend(loc='lower right')
        plt.plot([0, 1], [0, 1], 'r--')
        plt.ylabel('True Positive Rate')
        plt.xlabel('False Positive Rate')
        # plt.savefig('/Users/rohanpaul/Desktop/OnevAll.png'.format(cell), bbox_inches='tight')
        plt.show()
        plt.close()


singleCellLine(binary=True)
# OneVsAll()
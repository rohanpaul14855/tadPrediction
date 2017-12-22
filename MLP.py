import tensorflow as tf
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os, csv


def readData(cell):
    dirpath = '/Users/rohanpaul/Dropbox/AyLab/finalFeaturefilesLocal/train/{}'.format(cell)
    fname = os.path.join(dirpath, '{}_pos_features.txt'.format(cell))
    postrain = pd.read_csv(fname, index_col=0, header=0, sep='\t')

    fname = os.path.join(dirpath, '{}_neg_features.txt'.format(cell))
    negtrain = pd.read_csv(fname, index_col=0, header=0, sep='\t')

    dirpath = '/Users/rohanpaul/Dropbox/AyLab/finalFeaturefilesLocal/test/{}'.format(cell)
    fname = os.path.join(dirpath, '{}_pos_testFeatures.txt'.format(cell))
    postest = pd.read_csv(fname, index_col=0, header=0, sep='\t')

    fname = os.path.join(dirpath, '{}_neg_testFeatures.txt'.format(cell))
    negtest = pd.read_csv(fname, index_col=0, header=0, sep='\t')

    return pd.concat([postrain, negtrain]), pd.concat([postest, negtest])


# Arguments: mnist dataset
# returns a smaller training dataset composed of the first
# 1000 images and labels from each class
def split_data(mnist):
    small_x, small_y = [], []
    label = [0] * 10
    for i in range(10):
        j, num = 0, 0
        label[i] = 1
        while num < 1000:
            if (mnist.train.labels[j] == label).all():
                small_x.append(mnist.train.images[j])
                small_y.append(mnist.train.labels[j])
                num += 1
            j += 1
        label = [0] * 10
    return small_x, small_y


# Arguments: 2d list accuracy, containing training and testing accuracy
# learning rate (alpha) and the number of iterations after which the accuracy
# is calculated
# plots the erorrs using matplotlib
def plot(accuracy, alpha, step):
    plt.plot(accuracy[:, 0], 'r', accuracy[:, 1], 'b')
    plt.xlabel('Iterations (x{})'.format(step))
    plt.ylabel('Accuracy')
    plt.title('Training and Testing Accuracy\nLearning Rate: {}'.format(alpha))
    plt.legend(['Testing Accuracy', 'Training Accuracy'], loc=0)
    plt.show()
    plt.close()


# Arguments: features (xin) and labels (yin) for which accuracy is calculated
# Returns: the accuracy with the model so far
def metrics(xin, yin):
    predicted = sess.run(output, feed_dict={x: xin})
    trueLabels = tf.equal(tf.argmax(predicted, 1), tf.argmax(y, 1))
    acc = tf.reduce_mean(tf.cast(trueLabels, tf.float32))
    result = sess.run(acc, feed_dict={x: xin, y: yin})
    return result


x = tf.placeholder(tf.float32, [None, 25])
y = tf.placeholder(tf.float32, [None, 1])

# input -> hidden layer
w1 = tf.Variable(tf.random_normal([25, 300]))
b1 = tf.Variable(tf.random_normal([1, 300]))
l1 = tf.nn.relu(tf.matmul(x, w1) + b1)

# hidden -> output layer
w2 = tf.Variable(tf.random_normal([300, 1]))
b2 = tf.Variable(tf.random_normal([1, 1]))
output = tf.nn.tanh(tf.matmul(l1, w2) + b2)

# Training
ce = tf.reduce_mean(-tf.reduce_sum(y * tf.log(output), reduction_indices=[1]))
optimizer = tf.train.AdamOptimizer(0.5)
optimizer = optimizer.minimize(ce)
sess = tf.Session()

accuracy = []
start = tf.global_variables_initializer()
sess.run(start)


# xs, ys = split_data(mnist)


def batch(train, size):
    num_pos = sum(train['Class'] == 1)
    num_neg = train.shape[0] - num_pos
    weights = [0.5 / num_pos for _ in range(num_pos)] + [0.5 / num_neg for _ in range(num_neg)]
    batch = train.sample(n=size, axis=0, weights=weights)
    xbatch = np.array(batch[train.columns.values[:-1]])
    ybatch = np.array(batch['Class']).reshape(size, 1)
    return xbatch, ybatch


for cell in ['Hela', 'HMEC', 'HUVEC', 'NHEK', 'K562']:
    train, test = readData(cell)

    Xtrain = train[train.columns.values[:-1]]
    ytrain = train['Class']
    Xtest = test[test.columns.values[:-1]]
    ytest = np.array(test['Class']).reshape(test['Class'].shape[0], 1)
    smallData = False  # set to true to use first 1000 images of each class
    ################using small dataset
    for epoch in range(100):
        xbatch, ybatch = batch(train, size=100)  # ALWAYS USE EVEN NUMBER
        sess.run(optimizer, feed_dict={x: xbatch, y: ybatch})
        if epoch % 50 == 0:
            e1 = metrics(Xtest, ytest)
            e2 = metrics(xbatch, ybatch)
            print(e1, e2)
            accuracy.append([e1, e2])
            # accuracy.append(e2)

# plt.plot(accuracy)
# plt.show()
plot(np.array(accuracy), alpha=0.7, step=50)

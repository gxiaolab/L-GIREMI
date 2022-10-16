import pandas as pd
import numpy as np
from sklearn import linear_model
from functools import partial


def ecdf(x):
    '''
    Calculate the empirical cumulative density function

    Parameters:
        x - a list-like object contains the values.

    Return: a cumulative density function based on the input data.
    '''
    x = np.array(x)
    x.sort()
    n = len(x)
    y = np.concatenate([[0], np.linspace(1/n, 1, n)])

    def childfunc(sample, x, y, sorted=True):
        if not sorted:
            asort = np.argsort(x)
            x = np.take(x, asort, 0)
            y = np.take(y, asort, 0)
        idx = np.searchsorted(x, sample)
        return y[idx]

    return partial(childfunc, x=x, y=y)


def glm_score(data, model = 'lm'):
    # data with one column label train
    data_dummy = pd.concat(
        [
            pd.get_dummies(
                data[['change_type', 'up_seq', 'down_seq']]
            ),
            data['allelic_ratio_diff']
        ],
        axis=1
    )
    train_data_dummy = data_dummy.loc[data.label.notna()]
    train_label = data.loc[data.label.notna(), 'label'].map(
        lambda x: 1 if x == 'train_edit' else 0
    )
    if model == 'lm':
        glm = linear_model.LinearRegression()
    elif model == 'logistic':
        glm = linear_model.LogisticRegression()
    else:
        raise ValueError('Wrong model: {}.'.format(model))

    glm.fit(train_data_dummy.values, train_label.values)

    # train_pred = glm.predict(train_data_dummy.values)

    # train_accuracy = (
    #     train_pred == train_label
    # ).sum() / len(train_label)

    if model ==  'lm':
        data_score = glm.predict(data_dummy.values)
    elif model == 'logistic':
        data_score = glm.predict_proba(data_dummy.values)[:, 1]
    else:
        raise ValueError('Wrong model: {}.'.format(model))

    result = data.assign(score=data_score)
    return result


def score_performance(tp, tn, score):
    tp.loc[:, 'id'] = tp.apply(
        lambda a: '{}:{}:{}:{}:{}'.format(
            a['type'], a['chromosome'], a['pos'],
            a['strand'], a['change_type']
        ),
        axis = 1
    )

    tn.loc[:, 'id'] = tn.apply(
        lambda a: '{}:{}:{}:{}:{}'.format(
            a['type'], a['chromosome'], a['pos'],
            a['strand'], a['change_type']
        ),
        axis = 1
    )

    score.loc[:, 'id'] = score.apply(
        lambda a: '{}:{}:{}:{}:{}'.format(
            a['type'], a['chromosome'], a['pos'],
            a['strand'], a['change_type']
        ),
        axis = 1
    )

    scoredf = pd.DataFrame(
        {
            'score': score['score'].sort_values(
                ascending = False
            ).unique()
        }
    )
    scoredf.loc[:, 'tp'] = scoredf['score'].apply(
        lambda a:
        len(np.intersect1d(
            score.loc[score['score'] >= a, 'id'].values,
            tp['id'].values
        ))
    )
    scoredf.loc[:, 'fp'] = scoredf['score'].apply(
        lambda a:
        len(np.intersect1d(
            score.loc[score['score'] >= a, 'id'].values,
            tn['id'].values
        ))
    )
    scoredf.loc[:, 'fn'] = tp.shape[0] - scoredf['tp']
    scoredf.loc[:, 'tn'] = tn.shape[0] - scoredf['fp']
    scoredf.loc[:, 'precision'] = scoredf['tp'] / (
        scoredf['tp'] + scoredf['fp']
    )

    scoredf.loc[:, 'recall'] = scoredf['tp'] / (
        scoredf['tp'] + scoredf['fn']
    )

    scoredf.loc[:, 'f1'] = 2 * scoredf['precision'] * scoredf['recall'] / (
        scoredf['precision'] + scoredf['recall']
    )

    scoredf.loc[:, 'sensitivity'] = scoredf['recall']
    scoredf.loc[:, 'specificity'] = scoredf['tn'] / (
        scoredf['tn'] + scoredf['fp']
    )

    scoredf.loc[:, 'tpr'] = scoredf['recall']
    scoredf.loc[:, 'fpr'] = scoredf['fp'] / (scoredf['fp'] + scoredf['tn'])
    scoredf
    scoredf.loc[:, 'max_f1'] = False
    scoredf.loc[scoredf['f1'].argmax(), 'max_f1'] = True
    return scoredf


########################################

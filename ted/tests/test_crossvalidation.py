#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Tue 13 Aug 2013
#   Initial build.
#

"""
Test ted.sdss.cutouts
"""

def test_max_indices():

    import numpy as np

    from ted.sdss.cutouts.crossvalidation import max_indices

    mat0 = np.array([
        [2, 1, 2, 0],
        [2, 0, 1, 1],
        [2, 1, 1, 2],
        [2, 1, 1, 2],
    ])

    # Define answers
    answer = (
        np.array([0, 0, 1, 2, 2, 3, 3]),
        np.array([0, 2, 0, 0, 3, 0, 3])
    )

    # Get the results
    result = max_indices(mat0)

    # print answer
    # print result

    # Test
    assert np.all(result[0] == answer[0])
    assert np.all(result[1] == answer[1])


def test_get_centroids():

    import numpy as np

    from ted.sdss.cutouts.crossvalidation import get_centroids

    # max locs:
    # [0, 0, 1, 2, 2, 3, 3]
    # [0, 2, 0, 0, 3, 0, 3]
    mat0 = [
        [2, 1, 2, 0],
        [2, 0, 1, 1],
        [2, 1, 1, 2],
        [2, 1, 1, 2],
    ]
    # max locs:
    # [0, 0, 1, 1, 2, 3, 3]
    # [0, 2, 0, 3, 0, 0, 3]
    mat1 = [
        [5, 4, 5, 3],
        [5, 4, 4, 5],
        [5, 3, 4, 4],
        [5, 4, 4, 5],
    ]
    # max locs:
    # [0, 0, 1, 1, 2, 3, 3]
    # [0, 2, 0, 3, 0, 0, 3]
    mat2 = [
        [8, 7, 8, 6],
        [8, 7, 7, 8],
        [8, 7, 6, 7],
        [8, 7, 7, 8],
    ]

    # Stack arrays depth wise
    cube = np.dstack((mat0, mat1, mat2))

    # Define answers
    # answer = [
    #     (
    #         np.array([0, 0, 1, 2, 2, 3, 3]),
    #         np.array([0, 2, 0, 0, 3, 0, 3])
    #     )
    # ] * 3
    answer = [[2.0, 1.0], [1.0, 1.0], [1.0, 1.0]]

    # Get the results
    result = get_centroids(cube)

    print answer
    print result

    # Test
    assert np.all(result == answer)


def test_confusion_matrix():

    import numpy as np

    from ted.sdss.cutouts.crossvalidation import confusion_cube

    N_samples = 1
    N_folds = 5

    labels = np.vstack([np.ones(N_samples).astype(bool) for i in range(N_folds)])
    classes = np.array([True, False])

    if 1:
        """ALWAYS True"""
        """Test is made to pass this case"""
        always = np.ones(N_samples).astype(bool)
        predictions = np.vstack([always for i in range(N_folds)])
    else:
        """NEVER True"""
        never = np.zeros(N_samples).astype(bool)
        predictions = np.vstack([never for i in range(N_folds)])

    # Testing internals of the test
    # assert labels.shape == (N_folds, N_samples)
    # assert predictions.shape == (N_folds, N_samples)

    coc = confusion_cube(predictions, labels, classes)
    assert coc.shape == (2, 2, N_folds)

    mean = coc.mean(axis=2)
    std = coc.std(axis=2)
    assert np.all(mean == np.array([[N_samples, 0], [0, 0]]))
    assert np.all(std == np.array([[0, 0], [0, 0]]))

    count_vector = coc.sum(axis=0).sum(axis=0)
    assert np.all(count_vector == np.ones(N_samples))

    acc_vector = (coc[0, 0, :] + coc[1, 1, :]) / count_vector
    assert np.all(acc_vector == np.ones(N_samples))

    acc_mean = acc_vector.mean()
    acc_std = acc_vector.std()
    assert acc_mean == 1.0
    assert acc_std == 0.0


if __name__ == '__main__':
    pass

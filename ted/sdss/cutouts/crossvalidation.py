#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Fri 28 Mar 2014
#   Initial build.
#

import os
import time

import numpy as np
import pandas as pd
import yaml

# ... : ted
# from ... import msg
from ... import env
# from ... import TEDError
from ... import Namespace

# . : ted.sdss.cutouts
from . import load_cutout_sequences


class CVHelper(object):
    """Loads parameters and makes various parameters and methods available"""

    def __init__(self):
        self._load_parameters()
        self._load_data()

    # Initialise
    # ----------

    def _load_parameters(self):
        self._params = load_parameters()

    def _load_data(self):
        self._css, self._targets = load_cutout_sequences()

    @property
    def _cv(self): return self._params.get('cv')

    @property
    def _cs(self): return self._params.get('cs')

    @property
    def _xp(self): return self._params.get('xp')

    # Experiment properties
    # ---------------------

    def set_exp(self, key=None):
        if self._xp.get(key) is None: return None
        getattr(self, '_set_exp_' + key)()

    def _set_exp_any(self):
        xpp = dict2numpy(self._xp.get('any'))
        xp = Namespace()
        xp.name = 'any'
        xp.sigmas = xpp.get('sigmas')
        xp.taus = xpp.get('taus')
        xp.N_sigmas = xp.sigmas.size
        xp.N_taus = xp.taus.size
        xp.gskw = dict(exp=xp.name, **xpp)
        self.xp = xp
        # self.test()

    # Fold management
    # ---------------

    @property
    def _folds(self):
        return get_folds(self._css, self._targets, self._cv.get('N_folds'))

    # Cross validation
    # ----------------

    # def cv(self):
    #     for n, data in enumerate(self._folds):
    #         self._set_fold_parameters(fold=n + 1, data=data)
    #         # Training
    #         self._cube_of_predictions()
    #         self._matrix_of_accuracies()
    #         self._max_entries(moa)
    #         # Testing
    #         self._vector_of_predictions()
    #         self._scalar_of_accuracy()
    #         # Save to disk
    #         self._save_fold_results()

    def cv(self):
        self._cv_time = time.time()
        getattr(self, '_cv_' + self.xp.name)()

    def _cv_any(self):
        for n, data in enumerate(self._folds):
            # Fold number
            fnum = n + 1
            # Unpack training and test fold (f: features; l: labels)
            (train_f, train_l), (test_f, test_l) = data
            try:
                # raise Exception

                # Training set
                cop = self._cop_any(train_f)
                moa = self._moa_any(cop, test_f)

                self._save_fold_results_any(moa=moa, ftype='train', fnum=fnum)

            except Exception as e:
                print '_cv_any:', fnum, 'train'
                print '_cv_any:', train_f.shape, train_l.shape
                print '_cv_any:', e.message

            try:
                # raise Exception

                # Test set
                cop = self._cop_any(test_f)
                moa = self._moa_any(cop, test_l)
                self._save_fold_results_any(moa=moa, ftype='test', fnum=fnum)

            except Exception as e:
                print '_cv_any:', fnum, 'test'
                print '_cv_any:', test_f.shape, test_l.shape
                print '_cv_any:', e.message

    # Grid search
    # -----------

    def _cube_of_predictions(self, features):
        """Store result of training with chosen experiment"""
        # self._cop = getattr(self, '_cop_' + self.xp.name)()
        return getattr(self, '_cop_' + self.xp.name)(features)

    def _cop_any(self, features):
        """Return result of training with experiment ANY"""
        N = features.shape[0]
        cop = np.zeros((self.xp.N_sigmas, self.xp.N_taus, N)).astype(bool)
        for k, cs in enumerate(features):
            # if cs.has_gs_prediction and cs.gs_prediction_time > self._cv_time:
            if cs.has_gs_prediction and cs.gs_prediction_time > self._cv_time - 2 * 3600.:
                cop[:, :, k] = cs.gs_prediction
            else:
                cs.load(**self._cs)
                cs.calculate_residuals()
                cop[:, :, k] = cs.gridsearch(**self.xp.gskw)
                cs.cleanup()
        return cop

    def _moa_any(self, cop, labels):
        """Store training accuracies for experiment ANY.
        A.k.a.: matrix of accuracies (moa)."""
        print '_moa_any:', 'CoP.shape:', cop.shape
        print '_moa_any:', 'Labels.shape:', labels.shape
        return (
            # Broadcast along the parameter axes
            (cop == labels[None, None, :])
            # Sum along the samples (CutoutSequences)
            .sum(axis=2).astype(float)
            # Divide by the number of training samples
            / labels.size)

    # Single experiment
    # -----------------

    def _vector_of_predictions(self, features):
        """Store result of testing with chosen experiment"""
        return getattr(self, '_vop_' + self.xp)(features)

    def _vop_any(self, features):
        """Return result of testing with experiment ANY"""
        vop = np.zeros(features.shape[0]).astype(bool)
        for i, cs in enumerate(features):
            cs.load(**self._cs)
            cs.calculate_residuals()
            vop[i] = cs.predict(**self._parkw)
            cs.cleanup()
        return vop

    def _soa_any(self, vop, labels):
        """Store test accuracy for experiment ANY.
        A.k.a.: scalar of accuracy (soa)."""
        return (vop == labels).sum().astype(float) / labels.size

    # Single experiment-prediction kwargs
    # -----------------------------------

    @property
    def _parkw(self):
        """Property that holds the kwargs for a single experiment prediction"""
        return getattr(self, '_parkw_' + self.xp.name)

    @property
    def _parkw_any(self, row_ix, col_ix):
        """Prediction kwargs for experiment ANY"""
        return dict(
            exp=self.xp.name,
            sigma=self.xp.sigmas[row_ix],
            tau=self.xp.taus[col_ix])

    # Fold I/O
    # --------

    @property
    def _opath(self):
        # Hard-coded path for now ...
        # Would be nice, if the cutout size
        # could be stored in parameters.yaml.
        return os.path.join(env.paths.get('cutouts'), '101x101', 'results')

    @property
    def _fn_fstr(self):
        # return 'moa_{ft}_E-{xp}_CV-{N}-{n}_BG-{bg}_Q-{qstr}.csv'
        return 'moa_{}_E-{}_CV-{}-{}_BG-{}_Q-{}.csv'

    def _fn_kw(self, ftype, fnum):
        """Return current state (self._fold)"""
        return (
            ftype,
            self.xp.name,
            self._cv.get('N_folds'),
            fnum,
            self._cs.get('bg_model'),
            # The quality used for now
            # is just the same as the one used for loading the cutouts.
            self._qstr(self._cs.get('quality')),
        )

    @staticmethod
    def _qstr(quality):
        return ''.join([str(q) for q in sorted(quality)])

    def _save_fold_results(self):
        getattr(self, '_save_fold_results_' + self.xp.name)()

    def _save_fold_results_any(self, moa, ftype, fnum):
        # fname = self._fn_fstr.format(**self._fn_kw(ftype=ftype, fnum=fnum))
        fname = self._fn_fstr.format(*self._fn_kw(ftype=ftype, fnum=fnum))
        ofname = os.path.join(self._opath, fname)
        print 'Saving fold results to:', ofname
        # print moa.shape
        df = pd.DataFrame(data=moa, index=self.xp.sigmas, columns=self.xp.taus)
        df.to_csv(ofname, index=True, header=True)
        return ofname

    # Results I/O
    # -----------

    def _filenames(self, ftype):
        import itertools as it
        N = self._cv.get('N_folds')
        args = (
            [ftype],
            [self.xp.name],
            [N],
            range(1, N + 1),
            ['median'],
            [self._qstr(self._cs.get('quality'))]
        )
        for t in it.product(*args, repeat=1):
            ifname = os.path.join(self._opath, self._fn_fstr.format(*t))
            print ifname
            yield ifname

    def _load_results(self):
        """
        Based on given parameters and chosen experiment, load all
        combinations of filenames possible with the given parameters.

        """
        getattr(self, '_load_results_' + self.xp.name)()

    def _load_results_any(self, ftype):
        """Load results from experiment ANY"""

        # Cube of accuracies
        N_folds = self._cv.get('N_folds')
        coa = np.zeros((self.xp.N_sigmas, self.xp.N_taus, N_folds))
        for i, ifname in enumerate(self._filenames(ftype)):
            df = pd.read_csv(ifname, index_col=[0])
            coa[:, :, i] = df.values
        return coa

    # Visualisation
    # -------------

    def plot(self):
        """Plot results for chosen experiment"""
        getattr(self, '_plot_results_' + self.xp.name)()

    def _plot_results_any(self):
        """Plot results for experiment ANY"""

        from scipy import stats
        import matplotlib as mpl
        # Select backend beforehand to make it run on the image servers
        mpl.use('pdf')
        import matplotlib.pyplot as plt

        from mplconf import mplrc
        # from mplconf import rmath

        # Load data

        coa_train = self._load_results_any(ftype='train')
        coa_test = self._load_results_any(ftype='test')

        # print coa_train.shape
        # raise SystemExit

        # Get CV information

        N_folds = self._cv.get('N_folds')

        # Extract parameter space

        """
        NOTE
        ----
        All though these values are extracted from the parameter file,
        the cross validation may not have been run for these values of
        sigma and tau.
            Alternatively, the code should take into account what range
        of hyper-parameters were used (saved along with the results).
        """

        sigmas = self.xp.sigmas
        taus = self.xp.taus

        N_sigmas = sigmas.size
        N_taus = taus.size

        # Prepare the data for plotting
        # -----------------------------

        # Plot of every moa for each training fold
        # Check: This is what is loaded directly from the files on disk.

        # On top of this plot, I need the locations of the best
        # accuracies all together.
        # Matrices of maximum accuracies (momas), where the True entries
        # are those, for which the corresponding location in the moa is max.

        # TRAIN
        momas_train = np.zeros((N_sigmas, N_taus, N_folds)).astype(bool)
        for i in range(N_folds):
            momas_train[:, :, i] = max_locs(coa_train[:, :, i])

        # TEST
        momas_test = np.zeros((N_sigmas, N_taus, N_folds)).astype(bool)
        for i in range(N_folds):
            momas_test[:, :, i] = max_locs(coa_test[:, :, i])

        # Calculate the locations of the assumed best accuracies
        # This is the centre-of-mass/centroids for each set of
        # entries with the maximum value.

        # List of centroids for the best accuracy

        # TRAIN
        centroids_train = [
            centroid(
                max_indices(coa_train[:, :, i])
            )
            for i in range(N_folds)
        ]

        # TEST
        centroids_test = [
            centroid(
                max_indices(coa_test[:, :, i])
            )
            for i in range(N_folds)
        ]

        # Use slices to create matrices with NaNs
        # everywhere except at the centroid index.

        # TRAIN
        coms_train = np.zeros((N_sigmas, N_taus, N_folds))
        for i, (rix, cix) in enumerate(centroids_train):
            coms_train[rix, cix, i] = 1

        # NaNs are not shown in imshow-plots.
        coms_train[coms_train == 0] = np.nan

        # TEST
        coms_test = np.zeros((N_sigmas, N_taus, N_folds))
        for i, (rix, cix) in enumerate(centroids_test):
            coms_test[rix, cix, i] = 1

        # NaNs are not shown in imshow-plots.
        coms_test[coms_test == 0] = np.nan

        if 0:
            # Density estimation of the highest-accuracy locations
            kernel_train = stats.gaussian_kde(coms_train.sum(axis=2))

            # Density estimation of the highest-acuracies for the sum of the test accuracies.
            kernel_test = stats.gaussian_kde(coa_test.sum(axis=2))

        # Plot settings
        # -------------

        mplrc('publish_digital')

        # Extent
        extent = [taus.min(), taus.max(), sigmas.min(), sigmas.max()]

        # Image settings for the accuracies
        imkw = dict(origin='lower', aspect='auto')
        imkw.update(extent=extent, interpolation='nearest')

        # Image settings for matrices of maximum accuracies
        momaskw = imkw.copy()
        momaskw.update(cmap=mpl.cm.binary_r, alpha=.8)

        # Image settings for centre-of-mass matrices
        comskw = imkw.copy()
        comskw.update(cmap=None)

        # Accuracies can only be between 0 and 100 percent
        imkw.update(vmin=0, vmax=1)

        # Colorbar settings
        cbkw = dict(extend='neither', drawedges=False)
        cbkw.update(orientation='horizontal')
        cbkw.update(pad=.1)

        # Plot accuracies for each training fold overplotted with the best-accuracy positions.
        # -----

        fkw = dict(sharex=True, sharey=True, figsize=(15, 6))
        fig, axes = plt.subplots(2, N_folds, **fkw)

        iter_axes = zip(
            axes.flat[:N_folds],
            axes.flat[N_folds:])

        for i, (ax_top, ax_bot) in enumerate(iter_axes):
            # im = ax_top.imshow(coa_train[:, :, i], **imkw)
            ax_top.imshow(coa_train[:, :, i], **imkw)

            # Show the locations of all entries having the maximum accuracy
            ax_bot.imshow(momas_train[:, :, i], **momaskw)

            # Plot a red color for the location of the centre-of-mass
            # of the maximum accuracy indices.
            com = np.zeros((N_sigmas, N_taus, 4))
            ix = (coms_train[:, :, i] == 1)
            com[ix, :] = np.array([1., .0, .0, 1.])
            com[~ix, :] = np.array([0., .0, .0, .0])
            ax_bot.imshow(com, **comskw)

            # Rest of the plot setup
            ax_bot.set_xlabel(r'$\tau$')
            ax_bot.set_xticks([.1, .5, 1.])
            if i == 0:
                ax_top.set_ylabel(r'$\sigma$')
                ax_bot.set_ylabel(r'$\sigma$')

        # Plot colorbar above the accuracies...

        # ax = fig.add_axes([.25, 1.2, .5, .0])
        # fig.colorbar(mappable=im, ax=ax, **cbkw)

        fig.tight_layout()
        qstr = self._qstr(self._cs.get('quality'))
        fname = 'moa_train_folds+best_Q-{}.pdf'.format(qstr)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)

        print 'CORRECT THE EXTENT !!!'

        # Plot accuracies for each test fold overplotted with the best-accuracy positions.
        # -----

        fkw = dict(sharex=True, sharey=True, figsize=(15, 6))
        fig, axes = plt.subplots(2, N_folds, **fkw)

        iter_axes = zip(
            axes.flat[:N_folds],
            axes.flat[N_folds:])

        for i, (ax_top, ax_bot) in enumerate(iter_axes):
            # im = ax_top.imshow(coa_train[:, :, i], **imkw)
            ax_top.imshow(coa_test[:, :, i], **imkw)

            # Show the locations of all entries having the maximum accuracy
            ax_bot.imshow(momas_train[:, :, i], **momaskw)

            # Plot a red color for the location of the centre-of-mass
            # of the maximum accuracy indices.
            com = np.zeros((N_sigmas, N_taus, 4))
            ix = (coms_test[:, :, i] == 1)
            com[ix, :] = np.array([1., .0, .0, 1.])
            com[~ix, :] = np.array([0., .0, .0, .0])
            ax_bot.imshow(com, **comskw)

            # Rest of the plot setup
            ax_bot.set_xlabel(r'$\tau$')
            ax_bot.set_xticks([.1, .5, 1.])
            if i == 0:
                ax_top.set_ylabel(r'$\sigma$')
                ax_bot.set_ylabel(r'$\sigma$')

        # Plot colorbar above the accuracies...

        # ax = fig.add_axes([.25, 1.2, .5, .0])
        # fig.colorbar(mappable=im, ax=ax, **cbkw)

        fig.tight_layout()
        fname = 'moa_test_folds+best_Q-{}.pdf'.format(qstr)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)

        print 'CORRECT THE EXTENT !!!'

        raise SystemExit

        # Plot distribution of the locations of the best accuracies in the training folds.
        # -----
        fkw = dict(sharex=True, sharey=True, figsize=(15, 8))
        fig, (ax1, ax2) = plt.subplots(1, 2, **fkw)

        X, Y = np.mgrid[0:1:N_sigmas * 1j, 0:1:N_taus * 1j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        Z = np.reshape(kernel_train(positions).T, X.shape)
        ax1.imshow(np.rot90(Z), **imkw)

        # Rest of the plot setup
        ax1.set_ylabel(r'$\sigma$')
        for ax in (ax1, ax2):
            ax.set_xlabel(r'$\tau$')

        # Plot colorbar above the accuracies...
        # ax = fig.add_axes([.25, 1.2, .5, .0])
        # fig.colorbar(mappable=im, ax=ax, **cbkw)

        print 'HI'

        fig.tight_layout()
        fname = 'kde_coms+coa_test.sum.pdf'
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)

        print 'CORRECT THE EXTENT !!!'

        # Plot distribution of the sum of the test-fold accuracy matrices.
        # -----

        # Plot the surface of average accuracy for the test folds along with standard-deviation surfaces.
        # -----



    # -- END class CVHelper --


def max_indices(arr):
    """Get indices of the maximum"""
    return (arr == arr.max()).nonzero()


def max_locs(arr):
    """
    Get boolean array where True means
    entry value was same as maximum value.
    """
    return (arr == arr.flatten().max())


def centroid(indices_list):
    print 'Indices:'
    centroid = []
    for indices in indices_list:
        print indices
        centroid.append(np.round(np.mean(indices)))
    print ''
    print 'Centroid:', centroid
    print ''

    return centroid


def load_parameters(ifname=env.files.get('params')):
    """Load parameters for loading the data and running an experiment."""
    with open(ifname, 'r') as fsock:
        doc = yaml.load(fsock.read())
    return doc


def dict2numpy(d):
    """
    Walk through given dictionary recursively.

    If any sub-dictionary has keys that correspond to argument names
    for NunmPy functions such as

        numpy.linspace(start, stop, [num])  # here: start, stop, num
        numpy.arange([start,] stop[, step]) # here: start, stop, step

    the dictionary is replaced by the corresponding construct, i.e. either of

        numpy.linspace(**d)
        numpy.arange(**d)

    All three arguments must be specified in the dictionary for the parser to
    manipulate them.

    Parameters
    ----------
    d : dict
        A dictionary whose keys are to be parsed recursively.

    Returns
    -------
    d : dict
        A dictionary with any sub-dictionaries parsed recursively.
    L : numpy.linspace(**d)
    A : numpy.arange(**d)

    """

    if 3 == len(d.keys()):

        keys = ('start', 'stop',)

        if 3 == sum([d.has_key(k) for k in keys + ('num',)]):
            return np.linspace(**d)

        elif 3 == sum([d.has_key(k) for k in keys + ('step',)]):
            return np.arange(**d)

    for key, val in d.items():
        if isinstance(val, dict):
            d[key] = dict2numpy(val)

    return d


def get_folds(features, labels, N_folds):
    """Generator for easy loop over folds in cross validation."""

    features = np.array(features)
    labels = np.array(labels)

    # This seems reasonable to assume
    assert features.shape[0] == labels.shape[0]

    # For this project, I make sure that the
    # number of entries is divisible by `N_folds`
    assert features.shape[0] % N_folds == 0

    # Define chunk sizes
    # ------------------

    # Create training and test data sets
    N_features = features.shape[0]
    # Number of items in the test set
    # N_items is the number of samples per chunk
    N_items = N_features / N_folds

    # Define slices that will cut out the parts for training and testing
    # -------------
    tests = [
        slice(i * N_items, (i + 1) * N_items, None)
        for i in range(N_folds)]
    trains = [(
        slice(0, i * N_items, None),
        slice((i + 1) * N_items, N_folds * N_items, None))
        for i in range(N_folds)]

    # Generate the iterator
    for (s, (s1, s2)) in zip(tests, trains):

        test = (features[s], labels[s])

        train_features = np.append(features[s1], features[s2])
        train_labels = np.append(labels[s1], labels[s2])
        train = (train_features, train_labels)

        yield train, test


def LoG_radius2sigma(radius):
    if isinstance(radius, (int, float)):
        return radius / np.sqrt(2.)
    elif isinstance(radius, (list, np.ndarray)):
        return np.array(radius) / np.sqrt(2.)
    else:
        raise TypeError('`radius` is neither array-like or a number.')


def cv(exp='any'):  # , quality_combinations=[]):
    """Perform N-fold cross-validated grid search for chosen experiment."""
    cvh = CVHelper()
    cvh.set_exp('any')
    cvh.cv()


def analyse(exp='any'):  # , quality_combinations=[]):
    """Analyse results of N-fold cross-validated grid search for chosen experiment."""
    cvh = CVHelper()
    cvh.set_exp('any')
    cvh.analyse()
    # cvh.plot()


def plot(exp='any'):  # , quality_combinations=[]):
    """Plot results of N-fold cross-validated grid search for chosen experiment."""
    cvh = CVHelper()
    cvh.set_exp('any')
    cvh.plot()


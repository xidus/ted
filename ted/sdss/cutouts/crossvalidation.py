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
        self.N_folds = self._cv.get('N_folds')

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
    #     getattr(self, '_set_exp_' + key)()

    # def _set_exp_any(self):
    #     xpp = dict2numpy(self._xp.get('any'))
        xpp = dict2numpy(self._xp.get(key))
        xp = Namespace()
        xp.name = key
        xp.sigmas = xpp.get('sigmas')
        xp.taus = xpp.get('taus')
        xp.N_sigmas = xp.sigmas.size
        xp.N_taus = xp.taus.size
        xp.gskw = dict(exp=xp.name, **xpp)
        self.xp = xp

    @property
    def quality(self):
        if not hasattr(self, '_quality'):
            return self._cs.get('quality')
        else:
            return self._quality

    def set_quality(self, quality):
        if isinstance(quality, list):
            self._quality = quality

    @staticmethod
    def _qstr(quality):
        return ''.join([str(q) for q in sorted(quality)])

    @property
    def qstr(self):
        return self._qstr(self.quality)

    @property
    def _fname_prediction(self):
        return 'prediction_Q-{}.csv'.format(self.qstr)

    # Fold management
    # ---------------

    @property
    def _folds(self):
        return get_folds(self._css, self._targets, self.N_folds)

    # Cross validation
    # ----------------

    def cv(self):
        self._cv_time = time.time()
        getattr(self, '_cv_' + self.xp.name)()

    # ANY

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
                moa = self._moa(cop, train_l)
                self._save_fold_results_any(moa=moa, ftype='train', fnum=fnum)

            except Exception as e:
                print '_cv_any:', fnum, 'train'
                print '_cv_any:', e.message
                print '_cv_any:', '-' * 50, '\n'

            try:
                # raise Exception

                # Test set
                cop = self._cop_any(test_f)
                moa = self._moa(cop, test_l)
                self._save_fold_results_any(moa=moa, ftype='test', fnum=fnum)

            except Exception as e:
                print '_cv_any:', fnum, 'test'
                print '_cv_any:', e.message
                print '_cv_any:', '-' * 50, '\n'

    # BASELINE

    def _cv_blr(self): self._cv_baseline()
    def _cv_bla(self): self._cv_baseline()
    def _cv_bln(self): self._cv_baseline()

    def _cv_baseline(self):
        for n, data in enumerate(self._folds):
            # Fold number
            fnum = n + 1
            # Unpack training and test fold (f: features; l: labels)
            (train_f, train_l), (test_f, test_l) = data

            # Training set
            cop = self._cop_blr(train_f)
            moa = self._moa(cop, train_l)
            self._save_fold_results_any(moa=moa, ftype='train', fnum=fnum)

            # Test set
            cop = self._cop_blr(test_f)
            moa = self._moa(cop, test_l)
            self._save_fold_results_any(moa=moa, ftype='test', fnum=fnum)

    # Grid search
    # -----------

    # Cube of Predictions (CoP)

    # def _cube_of_predictions(self, features):
    #     """Store result of training with chosen experiment"""
    #     # self._cop = getattr(self, '_cop_' + self.xp.name)()
    #     return getattr(self, '_cop_' + self.xp.name)(features)

    # ANY

    def _cop_any(self, features):
        """Return result of training with experiment ANY"""
        shape = (self.xp.N_sigmas, self.xp.N_taus, features.shape[0])
        cop = np.zeros(shape).astype(bool)
        for k, cs in enumerate(features):
            cs.set_fname_gsp(self._fname_prediction)
            if cs.has_gs_prediction and cs.gs_prediction_time > self._cv_time:
            # if cs.has_gs_prediction and cs.gs_prediction_time > self._cv_time - 2 * 3600.:
                cop[:, :, k] = cs.gs_prediction
            else:
                cs.load(**self._cs)
                # print 'cs.quality:', cs.quality
                cs.set_quality(self.quality)
                # print 'cs.quality:', cs.quality
                cs.calibrate()
                # print 'cs.quality:', cs.quality
                # raise SystemExit
                cs.calculate_residuals()
                cop[:, :, k] = cs.gridsearch(**self.xp.gskw)
                cs.cleanup()
        return cop

    # BASELINE

    def _cop_blr(self, features):
        """Return CoP for baseline-experiment RANDOM (BLR)"""
        return self._cop_baseline(features)

    def _cop_bla(self, features):
        """Return CoP for baseline-experiment ALWAYS (BLA)"""
        return self._cop_baseline(features)

    def _cop_bln(self, features):
        """Return CoP for baseline-experiment NEVER (BLN)"""
        return self._cop_baseline(features)

    def _cop_baseline(self, features):
        """Return CoP for baseline experiment"""
        shape = (self.xp.N_sigmas, self.xp.N_taus, features.shape[0])
        cop = np.zeros(shape).astype(bool)
        for k, cs in enumerate(features):
            cs.set_fname_gsp(self._fname_prediction)
            if cs.has_gs_prediction and cs.gs_prediction_time > self._cv_time:
            # if cs.has_gs_prediction and cs.gs_prediction_time > self._cv_time - 2 * 3600.:
                cop[:, :, k] = cs.gs_prediction
            else:
                cop[:, :, k] = cs.gridsearch(**self.xp.gskw)
        return cop

    # Matrix of accuracies

    def _moa(self, cop, labels):
        """Store training accuracies for experiment ANY.
        A.k.a.: matrix of accuracies (moa)."""
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
            cs.set_quality(self.quality)
            cs.calibrate()
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
        opath = os.path.join(env.paths.get('cutouts'), '101x101', 'results')
        if not os.path.exists(opath):
            os.makedirs(opath)
        return opath

    @property
    def _fn_fstr(self):
        # return 'moa_{ft}_E-{xp}_CV-{N}-{n}_BG-{bg}_Q-{qstr}.csv'
        return 'moa_E-{}_Q-{}_{}_CV-{}-{}_BG-{}.csv'

    def _fn_kw(self, ftype, fnum):
        """Return current state (self._fold)"""
        return (
            self.xp.name,
            self._qstr(self.quality),
            ftype,
            self.N_folds,
            fnum,
            self._cs.get('bg_model'),
        )

    def _save_fold_results(self):
        getattr(self, '_save_fold_results_' + self.xp.name)()

    def _save_fold_results_any(self, moa, ftype, fnum):
        # fname = self._fn_fstr.format(**self._fn_kw(ftype=ftype, fnum=fnum))
        fname = self._fn_fstr.format(*self._fn_kw(ftype=ftype, fnum=fnum))
        ofname = os.path.join(self._opath, fname)
        # print 'Saving fold results to:', ofname
        print 'Saving fold results to:', fname
        # print moa.shape
        df = pd.DataFrame(data=moa, index=self.xp.sigmas, columns=self.xp.taus)
        df.to_csv(ofname, index=True, header=True)
        return ofname

    # Results I/O
    # -----------

    def _filenames(self, ftype):
        import itertools as it
        args = (
            [self.xp.name],
            [self._qstr(self.quality)],
            [ftype],
            [self.N_folds],
            range(1, self.N_folds + 1),
            [self._cs.get('bg_model')],
        )
        for t in it.product(*args, repeat=1):
            ifname = os.path.join(self._opath, self._fn_fstr.format(*t))
            print ifname
            yield ifname

    def _load_results(self, ftype):
        """
        Based on given parameters and chosen experiment, load all
        combinations of filenames possible with the given parameters.

        """
        #     getattr(self, '_load_results_' + self.xp.name)()

        # def _load_results_any(self, ftype):
        #     """Load results from experiment ANY"""

        # Cube of accuracies
        N_folds = self.N_folds
        coa = np.zeros((self.xp.N_sigmas, self.xp.N_taus, N_folds))
        for i, ifname in enumerate(self._filenames(ftype)):
            df = pd.read_csv(ifname, index_col=[0])
            coa[:, :, i] = df.values
        return coa

    def _load_prediction_cubes(self, ftype=None):
        """
        Return prediction cube for each fold and for each fold type.
        """

        if ftype is None: return None

        cops = []
        labels = []

        for n, data in enumerate(self._folds):

            # Unpack training and test fold (f: features; l: labels)
            (train_f, train_l), (test_f, test_l) = data

            if ftype == 'test':
                features = test_f
                labels.append(test_l)
            elif ftype == 'train':
                features = train_f
                labels.append(train_l)

            # Create new cube-of-predictions container
            shape = (self.xp.N_sigmas, self.xp.N_taus, features.shape[0])
            cop = np.zeros(shape).astype(bool)

            for k, cs in enumerate(features):
                cs.set_fname_gsp(self._fname_prediction)
                cop[:, :, k] = cs.gs_prediction

            cops.append(cop)

        return cops, labels

    # Visualisation
    # -------------

    def plot(self):
        """Plot results for chosen experiment"""
    #     getattr(self, '_plot_results_' + self.xp.name)()

    # def _plot_results_any(self):
    #     """Plot results for experiment ANY"""

        from scipy import stats
        import matplotlib as mpl
        # Select backend beforehand to make it run on the image servers
        mpl.use('pdf')
        import matplotlib.pyplot as plt

        from mplconf import mplrc
        from mplconf import rmath

        # Load data

        coa_train = self._load_results(ftype='train')
        coa_test = self._load_results(ftype='test')

        # print coa_train.shape
        # raise SystemExit

        # Get CV information
        N_folds = self.N_folds

        # Get experiment information
        qstr = self._qstr(self.quality)

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
        # Matrices of maximum accuracies (momas)

        """True => accuracy for parameter-pair location is max."""

        momas_train = get_matrices_of_maximum_accuracy(coa_train)
        momas_test = get_matrices_of_maximum_accuracy(coa_test)

        """
        Calculate the locations of the assumed best accuracies
        This is the centre-of-mass/centroid for each set of
        entries with the maximum value.
        """

        # List of centroids for the best accuracy

        centroids_train = get_centroids(coa_train)
        centroids_test = get_centroids(coa_test)

        # Use slices to create matrices with NaNs
        # everywhere except at the centroid index.

        shape = (N_sigmas, N_taus)
        coms_train = get_centroid_stack(centroids_train, shape, fill=np.nan)
        coms_test = get_centroid_stack(centroids_test, shape, fill=np.nan)

        if 0:
            # Density estimation of the highest-accuracy locations
            kernel_train = stats.gaussian_kde(coms_train.sum(axis=2))

            # Density estimation of the highest-acuracies for the sum of the test accuracies.
            kernel_test = stats.gaussian_kde(coa_test.sum(axis=2))

        # Plot settings
        # -------------

        # mplrc('publish_digital')
        mplrc('publish_printed')

        def imshow_com(img, ax=None):
            """Plot centre of mass image"""
            if ax is None: ax = plt.gca()
            ix = (img == 1)
            com = np.zeros((N_sigmas, N_taus, 4))
            com[ix, :] = np.array([1., .0, .0, 1.])
            com[~ix, :] = np.array([0., .0, .0, .0])
            ax.imshow(com, **comskw)

        # Extent
        # With the bottom is sigmas.min(), imshow must have origin='lower'
        tmin, tmax = taus.min(), taus.max()
        dtau = (taus[1] - taus[0]) / 2.
        smin, smax = sigmas.min(), sigmas.max()
        dsig = (sigmas[1] - sigmas[0]) / 2.
        extent = [(tmin - dtau), (tmax + dtau), (smin - dsig), (smax + dsig)]

        # GENERIC image settings
        imkw = dict(origin='lower', aspect='auto', interpolation='nearest')
        imkw.update(extent=extent)

        # Image settings for matrices of accuracies
        moaskw = imkw.copy()
        moaskw.update(cmap=mpl.cm.coolwarm)
        # Accuracies can only be between 0 and 100 percent
        moaskw.update(vmin=0, vmax=1)

        # Image settings for matrices of maximum-accuracy indices
        momaskw = imkw.copy()
        momaskw.update(cmap=mpl.cm.binary)

        # Image settings for matrices of centre-of-mass index
        comskw = imkw.copy()
        comskw.update(cmap=None)

        # Figure settings (for plt.subplots)
        ncols = 5
        nrows = N_folds / ncols
        figh = [3.3, 5.7]
        figsize = (15, figh[nrows - 1])
        fkw = dict(sharex=True, sharey=True, figsize=figsize)
        # fkw = dict(figsize=(15, 6.5))

        # Figure-adjustement settings
        figb = [.15, .15]
        adjustkw = dict(left=.03, bottom=figb[nrows - 1], top=.95, right=.95)
        adjustkw.update(wspace=.10, hspace=.10)

        # Colorbar settings
        cbkw = dict(extend='neither', drawedges=False)
        cbkw.update(ticks=np.linspace(.0, 1., 3))
        cbkw.update(orientation='vertical')
        # rect = [left, bottom, width, height]
        # crect = [0.96, 0.52, 0.01, 0.43]
        rect_b = adjustkw['bottom']
        rect_h = adjustkw['top'] - adjustkw['bottom']
        crect = [0.96, rect_b, 0.01, rect_h]

        # Accuracies
        # ----------

        """TRAIN: Plot accuracies for each training fold"""

        # Create the figure and axes
        fig, axes = plt.subplots(nrows, ncols, **fkw)

        # Plot data on eaxh axis
        for i, ax in enumerate(axes.flat):
            im = ax.imshow(coa_train[:, :, i], **moaskw)

        scl_kw = dict(l=r'$\sigma$', b=r'$\tau$', bticks=[.0, .5, 1.])
        set_common_labels(axes, ncols, **scl_kw)
        ax.set_xticks([.0, .5, 1.])

        # Plot colorbar to the right of the accuracies
        cax = plt.axes(crect)
        fig.colorbar(mappable=im, cax=cax, **cbkw)

        # fig.tight_layout()
        fig.subplots_adjust(**adjustkw)
        fig.suptitle(rmath('Train - Q = [{}]'.format(qstr)))
        fname = 'moa_E-{}_Q-{}_CV-{}_train_folds.pdf'.format(
            self.xp.name, qstr, N_folds)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)
        plt.close(fig)


        """TRAIN: Plot best-accuracy positions with centroids"""

        # Create the figure and axes
        fig, axes = plt.subplots(nrows, ncols, **fkw)

        # Plot data on eaxh axis
        for i, ax in enumerate(axes.flat):
            # Show the locations of all entries having the maximum accuracy
            ax.imshow(momas_train[:, :, i], **momaskw)

            # Plot a red color for the location of the centre-of-mass
            # of the maximum accuracy indices.
            imshow_com(img=coms_train[:, :, i], ax=ax)

        set_common_labels(axes, ncols, **scl_kw)
        ax.set_xticks([.0, .5, 1.])

        fig.tight_layout()
        fig.suptitle(rmath('Train - Q = [{}]'.format(qstr)))
        fname = 'moa_E-{}_Q-{}_CV-{}_train_best.pdf'.format(
            self.xp.name, qstr, N_folds)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)
        plt.close(fig)

        """TEST: Plot accuracies for each test fold"""

        # Create the figure and axes
        fig, axes = plt.subplots(nrows, ncols, **fkw)

        # Plot data on eaxh axis
        for i, ax in enumerate(axes.flat):
            im = ax.imshow(coa_test[:, :, i], **moaskw)

        # Set the x label on the bottom-most axes only
        for ax in axes.flat[-ncols:]:
            ax.set_xlabel(r'$\tau$')
            ax.set_xticks([.0, .5, 1.])

        # Set the y label on the left-most axes only
        for ax in axes.flat[::ncols]:
            ax.set_ylabel(r'$\sigma$')

        # Plot colorbar to the right of the accuracies
        cax = plt.axes(crect)
        fig.colorbar(mappable=im, cax=cax, **cbkw)

        # fig.tight_layout()
        fig.subplots_adjust(**adjustkw)
        fig.suptitle(rmath('Test - Q = [{}]'.format(qstr)))
        fname = 'moa_E-{}_Q-{}_CV-{}_test_folds.pdf'.format(
            self.xp.name, qstr, N_folds)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)
        plt.close(fig)

        """TEST: Plot best-accuracy positions with centroids"""

        # Create the figure and axes
        fig, axes = plt.subplots(nrows, ncols, **fkw)

        # Plot data on eaxh axis
        for i, ax in enumerate(axes.flat):
            # Show the locations of all entries having the maximum accuracy
            ax.imshow(momas_test[:, :, i], **momaskw)

            # Plot a red color for the location of the centre-of-mass
            # of the maximum accuracy indices.
            imshow_com(img=coms_test[:, :, i], ax=ax)

        set_common_labels(axes, ncols, **scl_kw)
        ax.set_xticks([.0, .5, 1.])

        fig.tight_layout()
        fig.suptitle(rmath('Train - Q = [{}]'.format(qstr)))
        fname = 'moa_E-{}_Q-{}_CV-{}_test_best.pdf'.format(
            self.xp.name, qstr, N_folds)
        ofname = os.path.join(self._opath, fname)
        plt.savefig(ofname)
        plt.close(fig)

        # raise SystemExit
        return

        # Plot distribution of the locations of the best accuracies in the training folds.
        # -----
        fkw = dict(sharex=True, sharey=True, figsize=(15, 8))
        fig, (ax1, ax2) = plt.subplots(1, 2, **fkw)

        X, Y = np.mgrid[0:1:N_sigmas * 1j, 0:1:N_taus * 1j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        Z = np.reshape(kernel_train(positions).T, X.shape)
        ax1.imshow(np.rot90(Z), **moaskw)

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

    # Analysis
    # --------

    def analyse(self):
        """Analyse results for chosen experiment"""
        getattr(self, '_analyse_' + self.xp.name)()

    def _analyse_any(self): self._analyse_simple()

    def _analyse_blr(self): self._analyse_simple() # baseline()
    def _analyse_bla(self): self._analyse_simple() # baseline()
    def _analyse_bln(self): self._analyse_simple() # baseline()

    def _analyse_baseline(self):
        """Analyse results for BASELINE experiment"""
        pass

    def _analyse_simple(self):
        """Analyse results any experiment"""

        # Load cubes of accuracies
        # (one matrix-of-accuracy for each fold and fold type)
        coa_train = self._load_results(ftype='train')

        # List of centroids for the best accuracy

        # Get location of the best parameters for each fold
        centroids = get_centroids(coa_train)

        # Use these indices to extract the predictions for each test fold
        # ultimately using these to get, not just accuracies, but the entire
        # table of confusion based on the confusion matrix (which for binary
        # classification is the same).
        # Get list of prediction cubes, one for each fold
        cops, labels = self._load_prediction_cubes(ftype='test')

        # print labels
        # raise SystemExit

        # Get a list of vectors containing the predictions
        predictions = []
        for fold_ix, (sigma_ix, tau_ix) in enumerate(centroids):
            predictions.append(cops[fold_ix][sigma_ix, tau_ix, :])

        # Get cube of confusion matrices
        classes = np.array([True, False])
        coc = confusion_cube(predictions, labels, classes)
        print 'coc.shape =', coc.shape
        print 'coc[:, :, 0]:'
        print coc[:, :, 0]

        print labels[0]
        print predictions[0]
        # raise SystemExit

        mean = coc.mean(axis=2)
        std = coc.std(axis=2)

        # count_vector contains total number of entries in each fold
        count_vector = coc.sum(axis=0).sum(axis=0).astype(float)

        # Divide by total number of entries for given fold to get accuracy
        #   Assuming here that the class of interest is the one indexed
        #   first along the rows and columns of the confusion matrix,
        #   i.e. the class SN.
        acc_vector = (coc[0, 0, :] + coc[1, 1, :]) / count_vector
        acc_mean = acc_vector.mean()
        acc_std = acc_vector.std()

        # Print all results
        print '\nClass order in confusion matrix:', list(classes)

        print '\nMean table of confusion:'
        print mean
        print '\nStandard deviation:'
        print std

        print '\nMean accuracy (TP + TN):', acc_mean,
        print '+-', acc_std / np.sqrt(coc.shape[2]),
        print '(std = {})'.format(acc_std)

    # -- END class CVHelper --


def max_indices(arr):
    """Get indices of the maximum"""
    return (arr == arr.max()).nonzero()


def max_locs(arr):
    """
    Get boolean array where True entries
    mark location of maximum value
    """
    return (arr == arr.max())


def get_matrices_of_maximum_accuracy(cube):
    """
    Get stack of boolean matrices where
    True marks location of maximum value
    for each layer in the stack.
    """
    # momas = np.zeros_like(cube).astype(bool)
    # for i in range(cube.shape[2]):
    #     momas[:, :, i] = max_locs(cube[:, :, i])
    # return momas
    return (cube == cube.max(axis=0).max(axis=0)[None, None, :])


def centroid(indices_list):
    """Return indices nearest to centre-of-mass of given indices"""
    return [np.round(np.mean(indices)) for indices in indices_list]


def get_centroids(cube):
    """Return list of centroids depth-wise for a data cube"""
    return [centroid(max_indices(cube[:, :, i]))
            for i in range(cube.shape[2])]


def get_centroid_stack(centroids, shape, fill=np.nan):
    """
    Create data cube with ones at the location
    of the centroids and `fill` elsewhere.
    """
    coms = np.zeros(shape + (len(centroids),))
    for i, (rix, cix) in enumerate(centroids):
        coms[rix, cix, i] = 1
    coms[coms == 0] = fill
    return coms


def set_common_labels(axes, ncols, l=None, b=None, lticks=None, bticks=None):

    if b is not None:
        # Set the x label on the bottom-most axes
        for ax in axes.flat[-ncols:]:
            ax.set_xlabel(b)

    if l is not None:
        # Set the y label on the left-most axes
        for ax in axes.flat[::ncols]:
            ax.set_ylabel(l)


def load_parameters(ifname=env.files.get('params')):
    """Load parameters for loading the data and running an experiment"""
    with open(ifname, 'r') as fsock:
        return yaml.load(fsock.read())


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


# -----------------------------------------------------------------------------

def table_of_confusion():
    """
                  Prediction
                  ----------

                  T         F

             ---------- ----------
    A|      |          |          |
    c|  T   |    TP    |    FN    |
    t|      |          | (Type 2) |
    u|       ---------- ----------
    a|      |          |          |
    l|  F   |    FP    |    TN    |
            | (Type 1) |          |
             ---------- ----------

    """
    pass


def _confusion_matrix(predictions, labels, classes=None):
    """

    Parameters
    ----------
    classes : 1D-array
        List of classes to compare against.
        If provided, the order of the entries is preserved along
        the rows and columns of the returned confusion matrix.

    Returns
    -------
    CFM : 2D-array
        The confusion matrix of shape (N, N) for N classes.
        The confusion matrix corresponding to every possible
        class that appears in `labels`.

    """

    import itertools as it

    predictions = np.asarray([predictions]).flatten()
    labels = np.asarray([labels]).flatten()
    if classes is None:
        # np.unique also sorts the data
        classes = np.unique(labels)
    else:
        pass
        # Can not use np.unique, since this also sorts the input array.
        # Can not use Python's set() function since its elements are unordered.
        # classes = np.array(list(set(classes)))
    N = classes.size

    # The order (0, 1) of the arguments for zip() translates into
    # (row-wise, col-wise) location of the two inputs. In this case,
    # since labels are given before predictions, the labels will be seen
    # along the rows, and the predictted classes shown  along the columns.
    z = zip(labels, predictions)
    p = it.product(classes, repeat=2)
    # print list(p)
    return np.array([z.count(x) for x in p]).reshape(N, N)


def confusion_cube(predictions, labels, classes=None):
    from sklearn.metrics import confusion_matrix as cm
    return np.dstack([cm(l, p, classes) for (p, l) in zip(predictions, labels)])


# -----------------------------------------------------------------------------


def cv(exp='any', quality=None):
    """Perform N-fold cross-validated grid search for chosen experiment."""
    cvh = CVHelper()
    cvh.set_exp(exp)
    if quality is not None:
        cvh.set_quality(quality)
    cvh.cv()


def analyse(exp='any', quality=None):
    """Analyse results of N-fold cross-validated grid search for chosen experiment."""
    cvh = CVHelper()
    cvh.set_exp(exp)
    if quality is not None:
        cvh.set_quality(quality)
    cvh.analyse()


def plot(exp='any', quality=None):
    """Plot results of N-fold cross-validated grid search for chosen experiment."""
    cvh = CVHelper()
    cvh.set_exp(exp)
    if quality is not None:
        cvh.set_quality(quality)
    cvh.plot()


# def analyse_baseline(quality=None):
#     """Analyse results of N-fold cross-validated grid search for BASELINE experiment."""
#     cvh = CVHelper()
#     if quality is not None:
#         cvh.set_quality(quality)

#     for exp in ('blr', 'bla', 'bln'):
#         cvh.set_exp(exp)
#         cvh.analyse()



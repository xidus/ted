#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: Thu 13 Mar 2014
#   Initial build.
#

"""
Plotting functions for CutoutSequence instances.

"""

import os

import numpy as np
import matplotlib as mpl
"""
Must declare the backend beforehand, since running the program on the
imageserver gives problems with the DISPLAY environment variable, when for some
reason saving as .pdf invokes the qt4 backend which requires a running
X server (which is not installed on the image(disk)servers) :/

This problem only occurs, when I run the program with Office Grid, not when I
generate the same figures using my IPython notebook server using the same code.

"""
mpl.use('pdf')
import matplotlib.pyplot as plt

# . : ted.sdss.cutouts
# .. : ted.sdss
from ..stripe82 import (
    ra_min, ra_max,
    # dec_min, dec_max,
    # stripe_width, stripe_height,
    stripe_extent_ra, stripe_extent_dec, w2h
)
# from . import load_cutout_sequences

from mplconf import mplrc
from mplconf import rmath

mplrc('publish_digital')


def get_extent(cs):
    """
    extent : None | (left, right, bottom, top)
    default : assigns zero-based row, column indices
              to the `x`, `y` centers of the pixels

    """
    return [0, cs.cube_shape[0] - 1, 0, cs.cube_shape[1] - 1]


def plot_covering(radec, fields, opath):
    """
    Visualise selection

    """

    # Get indices for those fields that have RA \in [300; 360]
    subtract_ix = fields.raMax > ra_max

    # Make copies of these arrays for the visualisation
    df_vis = fields.copy()
    radec_vis = radec.copy().flatten()

    # Translate the coordinate
    df_vis.raMax[subtract_ix] -= 360
    df_vis.raMin[subtract_ix] -= 360

    # Only if there wore any subtrations to be made,
    # subtract 360 from the chosen coordinate so that
    # its position in the stripe can be visualised.
    if subtract_ix.sum():
        radec_vis[0] -= 360

    # Prepare plot input
    pkw = dict(mec='None', mfc='w', ls='None', marker='o', ms=5)
    rkw = dict(fill=True, fc=mpl.rcParams['axes.color_cycle'][0], ec='None', alpha=.3)

    # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15., 15. * 8 / w2h * 2))
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15., 15. * 8 / w2h * 3))

    # TOP plot
    # For each covering field
    for i in range(df_vis.shape[0]):
        cov_ramin, cov_ramax = df_vis.raMin.iloc[i], df_vis.raMax.iloc[i]
        cov_decmin, cov_decmax = df_vis.decMin.iloc[i], df_vis.decMax.iloc[i]
        dra, ddec = (cov_ramax - cov_ramin), (cov_decmax - cov_decmin)
        ax1.add_patch(mpl.patches.Rectangle((cov_ramin, cov_decmin), dra, ddec, **rkw))
    ax1.plot([radec_vis[0]], [radec_vis[1]], **pkw)
    ax1.plot(stripe_extent_ra, stripe_extent_dec, c=mpl.rcParams['axes.color_cycle'][1])
    ax1.set_xlim(ra_max + 1, ra_min - 1)
    ax1.set_xlabel(r'$\alpha\ /\ \mathrm{deg}$')
    ax1.set_ylabel(r'$\delta\ /\ \mathrm{deg}$')

    # MIDDLE plot
    for i, (rmn, rmx) in enumerate(zip(df_vis.raMin, df_vis.raMax)):
        ax2.plot([rmn, rmx], [i] * 2)
    ax2.axvline(x=radec_vis[0], c='w')
    ax2.set_xlim(*ax2.get_xlim()[::-1])
    ax2.set_xlabel(r'$\alpha\ /\ \mathrm{deg}$')
    ax2.set_ylabel(rmath('Frame index'))
    # ax2.ticklabel_format(scilimits=(-10, 10))
    ax2.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))

    # BOTTOM plot
    for i, (dmn, dmx) in enumerate(zip(fields.decMin, fields.decMax)):
        ax3.plot([i] * 2, [dmn, dmx])
    ax3.axhline(y=radec_vis[1], c='w')
    ax3.set_xlabel(rmath('Frame index'))
    ax3.set_ylabel(r'$\delta\ /\ \mathrm{deg}$')

    fig.tight_layout()

    # Save it
    plt.savefig(os.path.join(opath, 'coverage.pdf'))
    # plt.savefig(os.path.join(opath, 'coverage.png'), dpi=72)

    plt.close(fig)


def plot_possible_cutouts(pxmax, opath):
    """
    POSIBLE Cutouts
        # cutout_overview()
    """

    # How many cutouts of a given pixel side length can be made
    # out of the covering frames of the given coordinate?

    # NB: These are not all the covering images from the initial selection (Selection I)
    #     They are only those images that allowed for the given cutout size, e.g. 101x101.

    nbins = 100

    # for i, pxsz in enumerate(pxmax):
    #     print i, pxsz

    fig, ax = plt.subplots(1, 1, figsize=(15., 4))

    h = ax.hist(pxmax, bins=np.linspace(0, 1489, nbins + 1), zorder=10)

    Hcum = np.cumsum(h[0][::-1])[::-1]
    Hoff = h[1][:-1]
    Hkw = dict(alpha=.5, fc=mpl.rcParams['axes.color_cycle'][1])
    ax.bar(left=Hoff, height=Hcum, width=Hoff[1] - Hoff[0], **Hkw)

    ax.set_ylabel(rmath('No. of frames with cutout size possible'))
    ax.set_xlabel(rmath('Cutout side length in pixels'))

    # ax.set_xlim(.0, 1500.)
    ax.set_xlim(ax.get_xlim()[0], 1500.)

    fig.tight_layout()

    plt.savefig(os.path.join(opath, 'max_cutout_size.pdf'))
    # plt.savefig(os.path.join(opath, 'max_cutout_size.png'), dpi=72)

    plt.close(fig)


def plot_pixel_indices(cs, opath):

    pkw = dict(ls='none', marker='o', ms=12, mec='None', alpha=.8)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15., 15. * 1489. / 2048. / 2))

    ax1.plot(cs.col_ixes_all, cs.row_ixes_all, **pkw)
    ax1.plot(cs.col_ixes, cs.row_ixes, **pkw)
    ax1.set_xlim(0, 2048 - 1)
    ax1.set_ylim(1489 - 1, 0)

    ax2.plot(cs.col_ixes_all_r, cs.row_ixes_all_r, **pkw)

    ax1.set_ylabel(rmath('y px coordinate'))
    for ax in (ax1, ax2):
        ax.set_xlabel(rmath('x px coordinate'))

    fig.tight_layout()

    plt.savefig(os.path.join(opath, 'pixel_indices_joint.pdf'))
    # plt.savefig(
    #     os.path.join(cs.path('coord'), 'pixel_indices_joint.png'), dpi=72)

    plt.close(fig)


def plot_time_coverage(cutout_dates, opath):

    # 2014-03-20:
    # This should not break the cutout-data creation loop.
    # If the data does dot support the operations below, skip this output.
    try:
        # Get matplotlib dates for the timeline
        cutout_mdates = mpl.dates.date2num(sorted(cutout_dates))
        cutout_mdates_diff = cutout_mdates[1:] - cutout_mdates[:-1]
        # cmdiff_mean = cutout_mdates_diff.mean()

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15., 4))

        ax1.axhline(y=1., c=mpl.rcParams['axes.color_cycle'][0])
        vkw = dict(ymin=.45, ymax=.55, lw=1.)
        vkw.update(c=mpl.rcParams['axes.color_cycle'][1])
        for mdate in cutout_mdates:
            ax1.axvline(x=mdate, **vkw)
        ax1.set_xlim(cutout_mdates.min(), cutout_mdates.max())
        ax1.set_ylim(.0, 2.)
        ax1.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y'))
        ax1.xaxis.set_major_locator(mpl.dates.YearLocator())
        ax1.set_yticklabels([])
        ax1.set_xlabel(rmath('Time'))

        # bins = [.0001, .001, .01, .1, 1, 10, 100, 1000, 10000]
        bins = np.logspace(-4, 4, 9)
        ax2.hist(cutout_mdates_diff, bins=bins)
        ax2.set_xscale('log')
        ax2.grid(b=False, which='minor')
        ax2.set_xlabel(rmath('Number of days between observations'))
        # print xticklabels

        # Create more human-readable axis
        fig.set_size_inches(15., 4 * 1.5)
        # fig.subplots_adjust(bottom=0.25)
        ax3 = ax2.twiny()
        ax3.set_frame_on(True)
        ax3.patch.set_visible(False)
        for spine in ax3.spines.itervalues():
             spine.set_visible(False)
        ax3b = ax3.spines['bottom']
        ax3b.set_visible(True)
        ax3b.set_position(('axes', -0.35))
        ax3b.axis.set_label_position('bottom')
        ax3b.axis.set_ticks_position('bottom')

        ax3.grid(b=False, which='both')
        ax3.set_xscale('log')
        locs = [1. / (24. * 60.), 1. / 24., 1., 30., 365., 365. * 5]
        lbls = ['1 m', '1 h', '1 d', '1 m', '1 yr', '5 yr']
        ax3.xaxis.set_ticks(locs)
        ax3.xaxis.set_ticks([], minor=True)
        ax3.xaxis.set_ticklabels(lbls)
        ax3.set_xlim(bins.min(), bins.max())
        ax3.set_xlabel(rmath('Time between observations'))

        fig.tight_layout()
        plt.savefig(os.path.join(opath, 'stats.pdf'))
        # plt.savefig(os.path.join(opath, 'stats.png'), dpi=72)

        plt.close(fig)

    except:
        err_msg = 'plot_time_coverage: AN ERROR OCCURRED !'
        print err_msg
        print 'plot_time_coverage: Skipping output for this CS instance ...'
        print 'plot_time_coverage: Saving incident in errors.log'
        ofname = os.path.join(opath, 'errors.log')
        with open(ofname, 'a') as fsock:
            fsock.write('{}\n'.format(err_msg))

        # There are not going to be many of these
        # (the first occurred almost halfway through tlist),
        # but since an open figure slows things down a bit,
        # it should be closed if possible.
        if 'fig' in dir():
            try:
                plt.close(fig)
            except:
                pass


def plot_background_models(cs):

    # Set common vmin and vmax
    # NOTE: I am using Python's built-in functions max and min, since
    #       NumPy's versions return nan for both np.max and np.min if
    #       this value is present in the array, whereas max and min
    #       do what I want in this case: return the minimum and maximum
    #       of the non-nans in the array.
    #
    vmin = min([min(bg_model.flat) for bg_model in cs.bg_models.values()])
    vmax = max([max(bg_model.flat) for bg_model in cs.bg_models.values()])

    # Choose overall colour map
    cmap = mpl.cm.cubehelix

    imkw = dict(cmap=cmap, interpolation='nearest')
    imkw.update(aspect='equal')  # , origin='lower')
    imkw.update(extent=get_extent(cs))

    cbkw = dict(extend='neither', drawedges=False)
    cbkw.update(pad=.005, orientation='horizontal')

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(15., 5.))

    imkw.update(vmin=np.log10(vmin), vmax=np.log10(vmax))
    im1 = ax1.imshow(np.log10(cs.bg_models.get('median')), **imkw)
    im2 = ax2.imshow(np.log10(cs.bg_models.get('mean')), **imkw)
    ax1.set_title(rmath('Log10(Median)'))
    ax2.set_title(rmath('Log10(Mean)'))

    imkw.update(vmin=vmin, vmax=vmax)
    im3 = ax3.imshow(cs.bg_models.get('median'), **imkw)
    im4 = ax4.imshow(cs.bg_models.get('mean'), **imkw)
    ax3.set_title(rmath('Median'))
    ax4.set_title(rmath('Mean'))

    # print np.log10(cs.bg_models.get('median')).flatten()[:10]
    # print cs.bg_models.get('median').flatten()[:10]

    tick_locator = mpl.ticker.MaxNLocator(nbins=5)
    auto_locator = mpl.ticker.AutoLocator()
    for ax, im in zip((ax1, ax2, ax3, ax4), (im1, im2, im3, im4)):
        ax.set_axis_off()
        cbar = fig.colorbar(mappable=im, ax=ax, **cbkw)
        cbar.locator = tick_locator
        cbar.ax.xaxis.set_major_locator(auto_locator)
        cbar.update_ticks()

    # for ax, im in zip((ax1, ax2), (im1, im2)):
    #     # cbar = fig.colorbar(mappable=im, ax=ax, **cbkw)
    #     fig.colorbar(mappable=im, ax=ax, **cbkw)

    # for ax, im in zip((ax3, ax4), (im3, im4)):
    #     # cbar = fig.colorbar(mappable=im, ax=ax, **cbkw)
    #     fig.colorbar(mappable=im, ax=ax, ticks=cticks, **cbkw)

    fig.tight_layout()
    plt.savefig(os.path.join(cs.path('coord'), 'bg_models.pdf'))

    plt.close(fig)


def check_offset(cs, offset=0):
    # ASSUMPTION: there are enough (MIN_NUMBER_OF_CUTOUTS) cutouts in
    # the sequence that offset will not be changed by both these lines.
    offset = max(offset - 2, 0)
    offset = min(offset, cs.cube_remap.shape[2] - 4)
    return offset


def plot_residual_samples(cs, offset=0):
    """
    Show a sample of difference images

    Parameters
    ----------
    offset : int
        Offset index to start from
        Use it to see the change between SN present and not

    """

    # Handle input
    offset = check_offset(cs, offset=offset)
    fname = 'residuals_offset-{:d}.pdf'.format(offset)
    ofname_check_residual = os.path.join(cs.path('results'), fname)

    cube = cs.cube_residual

    vmin = cube[:, :, offset:offset + 4].min()
    vmax = cube[:, :, offset:offset + 4].max()
    # print vmin, vmax

    # cmap = mpl.cm.bone
    cmap = mpl.cm.cubehelix

    imkw = dict(cmap=cmap, interpolation='nearest')
    imkw.update(aspect='equal')  # , origin='lower')
    imkw.update(vmin=vmin, vmax=vmax)

    cbkw = dict(extend='neither', drawedges=False)
    cbkw.update(pad=.005, orientation='horizontal')

    fig, axes = plt.subplots(1, 4, figsize=(15., 5.))

    for i in range(4):
        index = offset + i
        im = axes[i].imshow(cube[:, :, index], **imkw)
        axes[i].set_title('cube[:, :, {}]'.format(index))
        axes[i].set_axis_off()
        # cbar = fig.colorbar(mappable=im, ax=axes[i], **cbkw)
        fig.colorbar(mappable=im, ax=axes[i], **cbkw)

    fig.tight_layout()
    plt.savefig(ofname_check_residual)
    plt.close(fig)


def plot_LoG_samples(cs, offset=0):
    """
    Plot a couple of LoG-filtered images to check if filter size is
    reasonable (clear spot when SN present)

    Offset to start from
    Use it to see the change between SN present and not

    """

    # Handle input
    offset = check_offset(cs, offset=offset)
    fname = 'LoG_samples_offset-{:d}.pdf'.format(offset)
    ofname = os.path.join(cs.path('results'), fname)

    # Also plot them as 3D surfaces to get a better feel for the toplogy
    # of the the LoG-filtered difference images
    from mpl_toolkits.mplot3d import axes3d

    cube = cs.cube_LoG

    vmin = cube[:, :, offset:offset + 4].min()
    vmax = cube[:, :, offset:offset + 4].max()
    print 'vmin = {:.2f}\nvmax = {:.2f}'.format(vmin, vmax)

    cmap = mpl.cm.cubehelix

    # Imshow kwargs
    imkw = dict(cmap=cmap, interpolation='nearest')
    imkw.update(aspect='equal')  # , origin='lower')
    imkw.update(vmin=vmin, vmax=vmax)

    # Axes3D kwargs
    axkw = dict(cmap=cmap, rstride=5, cstride=5)
    axkw.update(vmin=vmin, vmax=vmax)

    # Colour bar kwargs
    cbkw = dict(extend='neither', drawedges=False)
    cbkw.update(pad=.005, orientation='horizontal')

    # Coordinates
    y = np.arange(cs.cube_shape[0])
    x = np.arange(cs.cube_shape[1])
    X, Y = np.meshgrid(x, y)

    # Plot
    fig = plt.figure(figsize=(15., 8))
    # fig, axes = plt.subplots(1, 4, figsize=(15., 5.))

    for i in range(4):

        index = offset + i
        j = i + 4
        # print i, j

        axi = fig.add_subplot(241 + i)
        axj = fig.add_subplot(241 + j, projection='3d')

        # In the top row
        im = axi.imshow(cube[:, :, index], **imkw)
        axi.set_title('cube_LoG[:, :, {}]'.format(index))
        axi.set_axis_off()
        # cbar = fig.colorbar(mappable=im, ax=axi, **cbkw)
        fig.colorbar(mappable=im, ax=axi, **cbkw)

        # In the bottom row
        Z = cube[:, :, index]
        # axj.plot_wireframe(X, Y, Z, **axkw)
        # surf = axj.plot_surface(X, Y, Z,  **axkw)
        axj.plot_surface(X, Y, Z,  **axkw)
        # cbar = fig.colorbar(mappable=surf, ax=axes[i], **cbkw)
        axj.set_zlim(vmin, vmax)

    fig.tight_layout()
    plt.savefig(ofname)
    plt.close(fig)


def plot_unclipped_cutouts(cs):
    pass


def plot_binary_fields(cs, offset=0):
    """
    Plot and save image of each binary field.

    Dependencies
    ------------
    cs.cube_minima_locs
    cs.cube_threshold
    tau : float
        relative threshold value
    quality : int, list
        list of qualities

    Needs
    -----
    I_min was one instance, change to others.

    """

    fname = 'check_minima.pdf'
    ofname_check_minima = os.path.join(cs.path('results'), fname)

    cmap = mpl.cm.binary

    imkw = dict(cmap=cmap, interpolation='nearest')
    imkw.update(aspect='equal')  # , origin='lower')
    imkw.update(extent=get_extent(cs))

    pkw = dict(ls='None', marker='o', ms=12, mec='None', alpha=.5)
    pkw.update(c=mpl.rcParams.get('axes.color_cycle')[1])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11., 5.))

    ax1.imshow(I_min, **imkw)
    ax1.plot([40], [40], **pkw)
    ax1.set_title(rmath('Minima'))

    ax2.imshow(I_min_thres, **imkw)
    ax2.plot([40], [40], **pkw)
    ax2.set_title(rmath('Minima (intensity > threshold)'))

    fig.tight_layout()
    plt.savefig(ofname_check_minima)
    plt.close(fig)


def plot_histogram_cutouts(cs):
    """
    For every cutout in cutout sequence

    Needs
    -----

    I_min was one instance, change to others.

    """

    # Where are the True entries?
    y_ix, x_ix = I_min.nonzero()

    # The intensities of the LoG-minima
    cube_masked = cube_remap[y_ix, x_ix, cut_ix]

    # Indices of intensities that are larger than zero.
    cix = cube_masked > 0

    # Histogram of pixel values for the original remapped image
    # pkw.update(c=mpl.rcParams.get('axes.color_cycle')[1])
    fig, ax = plt.subplots(1, 1, figsize=(15., 4.))
    h = ax.hist(cube_masked[cix], bins=50)
    # ax.set_yscale('log')
    ax.set_xlabel(rmath('Intensity'))
    ax.set_ylabel(rmath('Count'))
    ax.set_title(rmath('Pixel intensities from remapped frame'))
    fig.tight_layout()

    ofname_check_intensities = os.path.join(
        cs.path('results'), 'check_intensities_remap.pdf')
    plt.savefig(ofname_check_intensities)

    plt.close(fig)


def plot_histogram_intensity(cs):
    """
    For finding a suitable threshold

    Needs
    -----

    I_min was one instance, change to others.

    """

    # Where are the True entries?
    y_ix, x_ix = I_min.nonzero()

    # The intensities of the LoG-minima
    cube_masked = cube_remap[y_ix, x_ix, cut_ix]

    # Indices of intensities that are larger than zero.
    cix = cube_masked > 0

    # Histogram of pixel values for the original remapped image
    # pkw.update(c=mpl.rcParams.get('axes.color_cycle')[1])
    fig, ax = plt.subplots(1, 1, figsize=(15., 4.))
    h = ax.hist(cube_masked[cix], bins=50)
    # ax.set_yscale('log')
    ax.set_xlabel(rmath('Intensity'))
    ax.set_ylabel(rmath('Count'))
    ax.set_title(rmath('Pixel intensities from remapped frame'))
    fig.tight_layout()

    ofname_check_intensities = os.path.join(
        cs.path('results'), 'check_intensities_remap.pdf')
    plt.savefig(ofname_check_intensities)

    plt.close(fig)


def plot_histogram_LoG(cs):
    """
    How are the LoG pixel values distributed?

    """

    LoG = cube_LoG[:, :, cut_ix]
    # cix = LoG < np.inf
    cix = LoG < -1

    # Histogram of pixel values for the original remapped image
    # pkw.update(c=mpl.rcParams.get('axes.color_cycle')[1])
    fig, ax = plt.subplots(1, 1, figsize=(15., 4.))
    # h = ax.hist(array1D_LoG, bins=50)
    # h = ax.hist(cube_LoG[:, :, cut_ix].flat, bins=50)
    h = ax.hist(LoG[cix], bins=50)
    # ax.set_yscale('log')
    ax.set_title(rmath('Pixel intensities of the LoG-filtered residual'))
    fig.tight_layout()

    ofname_check_intensities = os.path.join(
        cs.path('results'), 'check_intensities_LoG.pdf')
    plt.savefig(ofname_check_intensities)

    plt.close(fig)


def plot_histogram_ranked_LoG(cs):
    """
    How are the relative changes in pixel intensity of
    the ranked and sorted LoG values distributed?

    """

    LoGs = LoG.copy().flatten()
    # print LoGs.shape
    LoGs.sort()
    # print LoGs[0], LoGs[-1]

    # Relative differences
    N_rel = 10
    LoGr = LoGs[:N_rel] / LoGs[1:N_rel + 1] - 1.

    print '; '.join(['{:.5%}'.format(R) for R in LoGr])


    if 0:
        # Histogram of pixel values for the original remapped image
        # pkw.update(c=mpl.rcParams.get('axes.color_cycle')[1])
        fig, ax = plt.subplots(1, 1, figsize=(15., 4.))
        # h = ax.hist(array1D_LoG, bins=50)
        # h = ax.hist(cube_LoG[:, :, cut_ix].flat, bins=50)
        h = ax.hist(LoG[cix], bins=50)
        # ax.set_yscale('log')
        ax.set_title(rmath('Pixel intensities of the LoG-filtered residual'))
        fig.tight_layout()

        # print h[0]

        ofname_check_intensities = os.path.join(opath, 'check_intensities_LoG.pdf')
        plt.savefig(ofname_check_intensities)

        plt.close(fig)


# def illustrate_cutout_sequences():

#     # Load all cutout sequences in tlist.
#     css, targets = load_cutout_sequences()

#     for cs in css:
#         plot_background_models(cs)
#         plot_residual_samples(cs)
#         plot_LoG_samples(cs)
#         plot_binary_fields(cs)


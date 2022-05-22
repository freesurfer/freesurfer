"""
This module provides a function to create screenshots

"""

# -----------------------------------------------------------------------------

def createScreenshots(SUBJECT, SUBJECTS_DIR, OUTFILE, INTERACTIVE = True, LAYOUT = None, 
    BASE = ["default"], OVERLAY = ["default"], SURF = ["default"], SURFCOLOR = ["default"], 
    VIEWS = ["default"]
    ):

    """
    function createScreenshots()
    
    Requires FREESURFER_HOME environment variable

    BASE, VIEWS must be lists, can be ["default"]

    OVERLAY, SURF, SURFCOLOR can be lists or None, can be ["default"]

    """

    # -----------------------------------------------------------------------------
    # auxiliary functions

    def computeLayout(n):
        
        import numpy as np
        
        y = np.ceil(np.sqrt(n))
        
        x = y - np.divmod(y**2 - n, y)[0]
        
        return int(x), int(y)

    # -----------------------------------------------------------------------------
    # imports

    import os
    import matplotlib

    import numpy as np
    import pandas as pd
    import nibabel as nb

    if not INTERACTIVE:
        matplotlib.use('Agg')

    from matplotlib import pyplot as plt
    from qatoolspython.qatoolspythonUtils import levelsetsTria

    # -----------------------------------------------------------------------------
    # settings

    FIGSIZE = 8

    FIGDPI = 100

    ALPHA = 0.5

    # -----------------------------------------------------------------------------
    # import image data

    if BASE == ["default"]:
        norm = nb.load(os.path.join(SUBJECTS_DIR, SUBJECT, 'mri', 'norm.mgz'))
    else:
        norm = nb.load(BASE[0])

    if OVERLAY is None:
        aseg = None
    elif OVERLAY == ["default"]:
        aseg = nb.load(os.path.join(SUBJECTS_DIR, SUBJECT, 'mri', 'aseg.mgz'))
    else:
        aseg = nb.load(OVERLAY[0])

    # -----------------------------------------------------------------------------
    # import surface data

    if SURF == ["default"]:
        surflist = [
            os.path.join(SUBJECTS_DIR, SUBJECT, 'surf', 'lh.white'),
            os.path.join(SUBJECTS_DIR, SUBJECT, 'surf', 'rh.white'),
            os.path.join(SUBJECTS_DIR, SUBJECT, 'surf', 'lh.pial'),
            os.path.join(SUBJECTS_DIR, SUBJECT, 'surf', 'rh.pial'),
            ]
    else:
        surflist = SURF

    surf = list()

    if surflist is not None:
        for i in range(len(surflist)):
            surf.append(nb.freesurfer.io.read_geometry(surflist[i], read_metadata=True))

    if SURFCOLOR == ["default"] and SURF == ["default"]:
        surfcolor = ['yellow','yellow','red','red']
    elif SURFCOLOR == ["default"] and SURF != ["default"]:
        surfcolor = ['yellow'] * len(surf)
    else:
        surfcolor = SURFCOLOR

    # -----------------------------------------------------------------------------
    # import colortable, compute auxiliary variables, and transform to matplotlib
    # colortable

    lut = pd.read_csv(os.path.join(os.environ.get('FREESURFER_HOME'),'FreeSurferColorLUT.txt'),
                      sep=' ',
                      comment='#',
                      header=None,
                      skipinitialspace=True,
                      skip_blank_lines=True,
                      error_bad_lines=False,
                      warn_bad_lines=True
                      )

    lut = np.array(lut)

    lutEnum = dict(zip(lut[:, 0], range(len(lut[:, 0]))))

    lutTab = np.array(lut[:, (2, 3, 4, 5)] / 255, dtype="float32")
    lutTab[:, 3] = 1

    lutMap = matplotlib.colors.ListedColormap(lutTab)

    # -----------------------------------------------------------------------------
    # determine VIEWS

    if VIEWS == ["default"]:
        CutsRRAS = [('x',-10),('x',10),('y',0),('z',0)]
    else:
        CutsRRAS = VIEWS

    # -----------------------------------------------------------------------------
    # compute levelsets

    LVL = list()

    # will not run if surf is empty (intended)
    for s in range(len(surf)):

        sLVL = list()

        for i in range(len(CutsRRAS)):
            
            # determine dimension
            if CutsRRAS[i][0] == 'x':
                iDim = 0
            elif CutsRRAS[i][0] == 'y':
                iDim = 1
            elif CutsRRAS[i][0] == 'z':
                iDim = 2

            # compute levelsets: e.g., LVL[SURF][VIEWS][vLVL|tLVL|iLVL][0][elementDimnnnnn1][elementDim2]
            sLVL.append(levelsetsTria(surf[s][0], surf[s][1], surf[s][0][:, iDim], CutsRRAS[i][1]))

        LVL.append(sLVL)

    # -----------------------------------------------------------------------------
    # get data for norm

    normData = norm.get_data()

    normVals = normData

    # -----------------------------------------------------------------------------
    # get data for aseg and change to enumerated aseg so that it can be used as
    # index to lutMap

    if aseg is not None:

        asegData = aseg.get_data()

        asegUnique, asegIdx = np.unique(asegData, return_inverse=True)

        asegEnum = np.array([lutEnum[x] for x in asegUnique])

        asegVals = np.reshape(asegEnum[asegIdx], (aseg.shape))

    # -----------------------------------------------------------------------------
    # compile image data for plotting

    # note that we are using the TkReg RAS system, not RAS per se
    # this is because the surfaces are in TkReg RAS also

    # one considerable advantage of the (TkReg) RAS system is that it is solely
    # based on switching dimensions, not rotating at odd angles etc, so no need
    # for interpolating.

    # x_R y_R z_R c_R
    # x_A y_A z_A c_A
    # x_S y_S z_S c_S
    #   0   0   0   1

    m = norm.header.get_vox2ras_tkr()
    n = norm.header.get_data_shape()

    xyzIdx = np.array(np.meshgrid(np.linspace(0, n[0] - 1, n[0]), np.linspace(0, n[1] - 1, n[1]), np.linspace(0, n[2] - 1, n[2])))
    xyzIdxFlat3 = np.reshape(xyzIdx, (3, np.prod(n))).transpose()
    xyzIdxFlat4 = np.hstack((xyzIdxFlat3, np.ones((np.prod(n), 1))))
    rasIdxFlat3 = np.matmul(m, xyzIdxFlat4.transpose()).transpose()[:, 0:3]

    normValsRAS = list()
    if aseg is not None:
        asegValsRAS = list()

    for i in range(len(CutsRRAS)):
        
        # determine dimension
        if CutsRRAS[i][0] == 'x':
            iDim = 0
        elif CutsRRAS[i][0] == 'y':
            iDim = 1
        elif CutsRRAS[i][0] == 'z':
            iDim = 2

        #
        sel = np.array(xyzIdxFlat3[rasIdxFlat3[:, iDim] == CutsRRAS[i][1], :], dtype="int")
        normValsRAS.append(np.squeeze(normVals[np.min(sel[:, 0]):np.max(sel[:, 0]) + 1, np.min(sel[:, 1]):np.max(sel[:, 1]) + 1, np.min(sel[:, 2]):np.max(sel[:, 2]) + 1]))
        if aseg is not None:
            asegValsRAS.append(np.squeeze(asegVals[np.min(sel[:, 0]):np.max(sel[:, 0]) + 1, np.min(sel[:, 1]):np.max(sel[:, 1]) + 1, np.min(sel[:, 2]):np.max(sel[:, 2]) + 1]))

    # -----------------------------------------------------------------------------
    # plotting: create a new figure, plot into it, then close it so it never gets
    # displayed

    # turn interactive plotting off unless interactive
    if not INTERACTIVE:
        plt.ioff()

    # compute layout
    if LAYOUT is None:
        myLayout = computeLayout(len(CutsRRAS))
    else:
        myLayout = LAYOUT

    # create list of de-facto layouts
    myLayoutList = list()
    for axsx in range(myLayout[0]):
        for axsy in range(myLayout[1]):
            myLayoutList.append((axsx,axsy))

    # create subplots
    fig, axs = plt.subplots(myLayout[0],myLayout[1])
    axs = np.reshape(axs,myLayout)

    # adjust layout
    fig.set_size_inches([FIGSIZE*myLayout[1],FIGSIZE*myLayout[0]])
    fig.set_dpi(FIGDPI)
    fig.set_facecolor('black')
    fig.set_tight_layout({'pad': 0})
    fig.subplots_adjust(wspace=0)

    # -----------------------------------------------------------------------------
    # plot each panel

    for p in range(len(CutsRRAS)):

        axsx = myLayoutList[p][0]
        axsy = myLayoutList[p][1]

        # determine dimensions
        if CutsRRAS[p][0] == 'x':
            # x axis of the image should be towards anterior, y axis should be towards superior in RAS image
            dims = (1, 2)
            # determine extent
            extent = (rasIdxFlat3[0, dims[0]], rasIdxFlat3[-1, dims[0]], rasIdxFlat3[0, dims[1]], rasIdxFlat3[-1, dims[1]])
            # imshow puts the first dimension (rows) of the data on the y axis, and the second (columns) on the x axis
            cor = np.where(m[dims[0], 0:3])[0]
            axi = np.where(m[dims[1], 0:3])[0]
            #
            if axi < cor:
                axs[axsx,axsy].imshow(normValsRAS[p], cmap='gray', origin='lower', extent=extent)
                if aseg is not None:
                    axs[axsx,axsy].imshow(asegValsRAS[p] + 0.5, cmap=lutMap, origin='lower', extent=extent, vmin=0, vmax=len(lutTab), alpha=ALPHA)
            else:
                axs[axsx,axsy].imshow(normValsRAS[p].transpose(), cmap='gray', origin='lower', extent=extent)
                if aseg is not None:
                    axs[axsx,axsy].imshow(asegValsRAS[p].transpose() + 0.5, cmap=lutMap, origin='lower', extent=extent, vmin=0, vmax=len(lutTab), alpha=ALPHA)
        elif CutsRRAS[p][0] == 'y':
            #
            dims = (0, 2)
            # determine extent
            extent = (rasIdxFlat3[0, dims[0]], rasIdxFlat3[-1, dims[0]], rasIdxFlat3[0, dims[1]], rasIdxFlat3[-1, dims[1]])
            # imshow puts the first dimension (rows) of the data on the y axis, and the second (columns) on the x axis
            sag = np.where(m[dims[0], 0:3])[0]
            axi = np.where(m[dims[1], 0:3])[0]
            #
            if axi < sag:
                axs[axsx,axsy].imshow(normValsRAS[p].transpose(), cmap='gray', origin='lower', extent=extent)
                if aseg is not None:
                    axs[axsx,axsy].imshow(asegValsRAS[p].transpose() + 0.5, cmap=lutMap, origin='lower', extent=extent, vmin=0, vmax=len(lutTab), alpha=ALPHA)
            else:
                axs[axsx,axsy].imshow(normValsRAS[p].transpose(), cmap='gray', origin='lower', extent=extent)
                if aseg is not None:
                    axs[axsx,axsy].imshow(asegValsRAS[p].transpose() + 0.5, cmap=lutMap, origin='lower', extent=extent, vmin=0, vmax=len(lutTab), alpha=ALPHA)
        elif CutsRRAS[p][0] == 'z':
            #
            dims = (0, 1)
            # determine extent
            extent = (rasIdxFlat3[0, dims[0]], rasIdxFlat3[-1, dims[0]], rasIdxFlat3[0, dims[1]], rasIdxFlat3[-1, dims[1]])
            # imshow puts the first dimension (rows) of the data on the y axis, and the second (columns) on the x axis
            sag = np.where(m[dims[0], 0:3])[0]
            cor = np.where(m[dims[1], 0:3])[0]
            if axi < sag:
                axs[axsx,axsy].imshow(normValsRAS[p].transpose(), cmap='gray', origin='lower', extent=extent)
                if aseg is not None:
                    axs[axsx,axsy].imshow(asegValsRAS[p].transpose() + 0.5, cmap=lutMap, origin='lower', extent=extent, vmin=0, vmax=len(lutTab), alpha=ALPHA)
            else:
                axs[axsx,axsy].imshow(normValsRAS[p].transpose(), cmap='gray', origin='lower', extent=extent)
                if aseg is not None:
                    axs[axsx,axsy].imshow(asegValsRAS[p].transpose() + 0.5, cmap=lutMap, origin='lower', extent=extent, vmin=0, vmax=len(lutTab), alpha=ALPHA)

        # prepare plot
        if rasIdxFlat3[0, dims[0]] > rasIdxFlat3[-1, dims[0]]:
            axs[axsx,axsy].invert_xaxis()
        if rasIdxFlat3[0, dims[1]] > rasIdxFlat3[-1, dims[1]]:
            axs[axsx,axsy].invert_yaxis()

        axs[axsx,axsy].set_axis_off()
        axs[axsx,axsy].set_aspect('equal')

        # now plot
        for s in range(len(surf)):
            for i in range(len(LVL[s][p][1][0])):
                axs[axsx,axsy].plot(
                    (LVL[s][p][0][0][LVL[s][p][1][0][i][0] - 1][dims[0]], LVL[s][p][0][0][LVL[s][p][1][0][i][1] - 1][dims[0]]),
                    (LVL[s][p][0][0][LVL[s][p][1][0][i][0] - 1][dims[1]], LVL[s][p][0][0][LVL[s][p][1][0][i][1] - 1][dims[1]]),
                    color=surfcolor[s],linewidth=np.round(FIGSIZE/8))

    # -----------------------------------------------------------------------------
    # output

    if not INTERACTIVE:
        plt.savefig(OUTFILE, facecolor=fig.get_facecolor())
        plt.close(fig)

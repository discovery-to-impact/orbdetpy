import numpy as np
import scipy.io as sio
import scipy.linalg
import json


def rtsSmoother(xPre, xPost, PPre, PPost, filetype):
    N = np.size(xPost, 0)
    T = np.size(xPost, 1)
    time = np.linspace(1, T, T)
    alpha = 1
    k = 3 - N
    wi_c = weight(N, alpha, k, filetype)
    xS = np.zeros([N, T])
    PS = np.zeros([N, N, T])

    for kk in range(T, 0, -1):

        if kk == T:
            xS[:, [kk - 1]] = xPost[:, [kk - 1]]
            PS[:, :, [kk - 1]] = PPost[:, :, [kk - 1]]
        else:
            SPpre = sigmapoints(xPre[:, [kk]], PPre[:, :, kk], alpha, k)
            SPpost = sigmapoints(xPost[:, [kk - 1]], PPost[:, :, kk - 1], alpha, k)
            C = np.zeros(N)
            for ii in range(0, 2 * N, 1):
                C = C + wi_c * (SPpost[:, [ii]] - xPost[:, [kk - 1]]) * np.transpose(SPpre[:, [ii]] - xPre[:, [kk]])

            A = np.dot(C, np.linalg.inv(PPre[:, :, kk]))
            xS[:, [kk - 1]] = xPost[:, [kk - 1]] + np.dot(A, (xS[:, [kk]] - xPre[:, [kk]]))
            PS[:, :, kk - 1] = PPost[:, :, kk - 1] + np.dot(np.dot(A, (PS[:, :, kk] - PPre[:, :, kk])), np.transpose(A))

    return xS, PS, time


def dataConfigMat(filename):
    filetype = 'mat'
    content = sio.loadmat(filename)
    xPre = content['store_xBar']
    xPost = content['store_xHat']
    PPre = content['store_PBar']
    PPost = content['store_P']

    return xPre, xPost, PPre, PPost, filetype


def dataConfigJSON(filename):
    filetype = 'json'
    with open(filename) as f:
        jsonData = json.load(f)

    dataSize0 = np.size(jsonData["Estimation"])
    dataSize1 = np.size(jsonData["Estimation"][0]["xBar"])
    xPre = np.zeros([dataSize1, dataSize0])
    xPost = np.zeros([dataSize1, dataSize0])
    PPre = np.zeros([dataSize1, dataSize1, dataSize0])
    PPost = np.zeros([dataSize1, dataSize1, dataSize0])
    for ii in range(0, dataSize0, 1):
        xPre[:, [ii]] = np.reshape(np.asarray(jsonData["Estimation"][ii]["xBar"]), [dataSize1, 1])
        xPost[:, [ii]] = np.reshape(np.asarray(jsonData["Estimation"][ii]["EstimatedState"]), [dataSize1, 1])
        PPre[:, :, [ii]] = np.reshape(np.asarray(jsonData["Estimation"][ii]["PBar"]), [dataSize1, dataSize1, 1])
        PPost[:, :, [ii]] = np.reshape(np.asarray(jsonData["Estimation"][ii]["EstimatedCovariance"]),
                                       [dataSize1, dataSize1, 1])

    return xPre, xPost, PPre, PPost, filetype


def weight(n, alpha, k, filetype):
    lambda_ = alpha ** 2 * (n + k) - n
    if filetype == 'json':
        wi_c = 1 / (2 * n)  # for JSON data
    else:
        wi_c = 1 / (2 * (n + lambda_))  # for mat data

    return wi_c


def sigmapoints(xBar, PBar, alpha, k):
    n = np.size(xBar, 0)
    lambda_ = alpha ** 2 * (n + k) - n
    L = scipy.linalg.cholesky(PBar, lower=True)
    sigmaPoints = np.empty([n, 2 * n])
    sigmaPoints[:, [0]] = xBar
    sigmaPoints[:, 1:n + 1] = xBar + np.sqrt(n + lambda_) * L
    sigmaPoints[:, n:2 * n + 1] = xBar - np.sqrt(n + lambda_) * L
    return sigmaPoints

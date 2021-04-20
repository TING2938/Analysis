import cv2
import numpy as np 

def RGBtoBGR(img):
    b = img[:, :, 0].copy()
    g = img[:, :, 1].copy()
    r = img[:, :, 2].copy()

    img[:, :, 0] = r
    img[:, :, 1] = g
    img[:, :, 2] = b
    
def toGray(img):
    out = img[:, :, 0]*0.0722+img[:, :, 1]*0.7152+img[:, :, 2]*0.2126
    return out.astype(np.uint8)

def thresholding(img):
    img[img < 128] = 0
    img[img >= 128] = 255

def otsuMethod(gray):
    t = np.arange(1, 255, 1)
    sqSb = np.zeros(t.shape)
    for ndx, t0 in enumerate(t):
        w0 = len(gray[gray < t0]) / len(gray)
        w1 = 1 - w0
        M0 = gray[gray < t0].mean()
        M1 = gray[gray >= t0].mean()
        sqSb[ndx] = w0 * w1 * (M0 - M1)**2
    minarg = sqSb.argmin()
    tSel = t[minarg]
    gray[gray < tSel] = 0
    gray[gray >= tSel] = 255


img = cv2.imread("imori.jpg")
cv2.imshow("img", img)

out = toGray(img)
# thresholding(out)
otsuMethod(out)
cv2.imshow("out", out)

cv2.waitKey(0)

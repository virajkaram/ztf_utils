import os, re

filtdict = {"vv": "V", "bb": "B", "uu": "U", "w1": "UW1", "m2": "UM2", 
            "w2": "UW2", "wh": "W"}

dir = os.listdir(".")

os.chdir("uvot/image")

dir2 = os.listdir(".")
for image in dir2:

    if re.search("_sk.img", image):
        print(image)
        filt = filtdict[image[14:16]]
        os.system("uvotimsum %s %s.fits" % (image, filt))
        os.system("cp ../../*.reg .")
        os.system("uvotsource image=%s.fits srcreg=src.reg bkgreg=bkg.reg sigma=3.0 outfile=%s.out syserr=yes output=ALL apercorr=CURVEOFGROWTH > %s.dat" % (filt, filt, filt))

os.chdir("../../../")
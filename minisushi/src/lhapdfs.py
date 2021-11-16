import json
import lhapdf

from init import *
from utils import *
from config import *

class sushi_pdfs:
    _pdfnamein = einital().initial()['pdfnamein']
    _muf = einital().initial()['muFggh']
    _ppcoll = einital().initial()['ppcoll']
    """
    class sushi_pdfs(enitial):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cinit = super().initial()
        self.pdfnamein = self.cinit['pdfnamein']
        self.muf = self.cinit['muFggh']
    """
    def initialize_pdf(self):
        lhapdf.setVerbosity(0)
        pdfnamein = self._pdfnamein
        p = lhapdf.mkPDF(pdfnamein)
        return p
    
    def pdfs(self, x):
        p = self.initialize_pdf()
        ff = [p.xfxQ(i, x, self._muf) for i in range(-5,6)]
        return ff
    
    def pdfs_batch(self, x):
        p = self.initialize_pdf()
        x = np.array(x)
        if(x.shape == ()):
            x = np.array([x])
        
        ff = [[p.xfxQ(i, x, self._muf) for i in range(-5,6)] for x in x]
        return np.array(ff)

    def getpdfs(self, x):
        temp = self.pdfs(x)

        bot = temp[10]
        chm = temp[9]
        stra = temp[8]
        up = temp[7]
        dn = temp[6]
        glu = temp[5]
        dsea = temp[4]
        usea = temp[3]
        sbar = temp[2]
        cbar = temp[1]
        bbar = temp[0]

        return up, dn, usea, dsea, stra, sbar, chm, cbar, bot, bbar, glu
    
    def getpdfs_batch(self, x):
        temp = self.pdfs_batch(x)
        bot = temp[:,10]
        chm = temp[:,9]
        stra = temp[:,8]
        up = temp[:,7]
        dn = temp[:,6]
        glu = temp[:,5]
        dsea = temp[:,4]
        usea = temp[:,3]
        sbar = temp[:,2]
        cbar = temp[:,1]
        bbar = temp[:,0]

        return up, dn, usea, dsea, stra, sbar, chm, cbar, bot, bbar, glu
    
    
    def PDFgg(self, x1, x2):
        upv1, dnv1, usea1, dsea1, stra1, sbar1, chm1, cbar1, bot1, bbar1, glu1 = self.getpdfs(x1)
        upv2, dnv2, usea2, dsea2, stra2, sbar2, chm2, cbar2, bot2, bbar2, glu2 = self.getpdfs(x2)
        

        PDFgg = glu1*glu2
        return PDFgg

    def PDFgg_batch(self, x1, x2):
        upv1, dnv1, usea1, dsea1, stra1, sbar1, chm1, cbar1, bot1, bbar1, glu1 = self.getpdfs_batch(x1)
        upv2, dnv2, usea2, dsea2, stra2, sbar2, chm2, cbar2, bot2, bbar2, glu2 = self.getpdfs_batch(x2)
        

        PDFgg = glu1*glu2
        return PDFgg

    def PDFqg(self, x1, x2):
        upv1, dnv1, usea1, dsea1, stra1, sbar1, chm1, cbar1, bot1, bbar1, glu1 = self.getpdfs(x1)
        upv2, dnv2, usea2, dsea2, stra2, sbar2, chm2, cbar2, bot2, bbar2, glu2 = self.getpdfs(x2)

        q1 = upv1 + dnv1+usea1 + dsea1 + stra1 + sbar1 + chm1 + cbar1 + bot1 + bbar1
        PDFqg = q1*glu2
        return PDFqg

    def PDFqg_batch(self, x1, x2):
        upv1, dnv1, usea1, dsea1, stra1, sbar1, chm1, cbar1, bot1, bbar1, glu1 = self.getpdfs_batch(x1)
        upv2, dnv2, usea2, dsea2, stra2, sbar2, chm2, cbar2, bot2, bbar2, glu2 = self.getpdfs_batch(x2)

        q1 = upv1 + dnv1+usea1 + dsea1 + stra1 + sbar1 + chm1 + cbar1 + bot1 + bbar1
        PDFqg = q1*glu2
        return PDFqg

    def PDFqq(self, x1, x2):
        upv1, dnv1, usea1, dsea1, str1, sbar1, chm1, cbar1, bot1, bbar1, glu1 = self.getpdfs(x1)
        upv2, dnv2, usea2, dsea2, str2, sbar2, chm2, cbar2, bot2, bbar2, glu2 = self.getpdfs(x2)

        if(self._ppcoll):
            # ppbar
            PDFqq = upv1*usea2+upv2*usea1+dnv1*dsea2+dnv2*dsea1+str1*sbar2+str2*sbar1+chm1*cbar2+chm2*cbar1+bot1*bbar2+bot2*bbar1
        else:
            # pp
            PDFqq=upv1*upv2+usea1*usea2+dnv1*dnv2+dsea1*dsea2+ 2*(str1*str2+chm1*chm2+bot1*bot2)
        return PDFqq

    def PDFqq_batch(self, x1, x2):
        upv1, dnv1, usea1, dsea1, str1, sbar1, chm1, cbar1, bot1, bbar1, glu1 = self.getpdfs_batch(x1)
        upv2, dnv2, usea2, dsea2, str2, sbar2, chm2, cbar2, bot2, bbar2, glu2 = self.getpdfs_batch(x2)

        if(self._ppcoll):
            # ppbar
            PDFqq = upv1*usea2+upv2*usea1+dnv1*dsea2+dnv2*dsea1+str1*sbar2+str2*sbar1+chm1*cbar2+chm2*cbar1+bot1*bbar2+bot2*bbar1
        else:
            # pp
            PDFqq=upv1*upv2+usea1*usea2+dnv1*dnv2+dsea1*dsea2+ 2*(str1*str2+chm1*chm2+bot1*bot2)
        return PDFqq

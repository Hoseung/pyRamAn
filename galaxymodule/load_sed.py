class Sed():
    def __init__(self, sed_blocks=[]):
        sb = sed_blocks[0]
        self.nmetals=len(sed_blocks)
        self.n_ages = sb.n_ages
        self.n_wavs = sb.n_wavs
        self.age_points = sb.age_points * 1e-9 # in Gyr
        self.sed_wavelengths = sb.wavelengths
        self.isochrone = sb.isochrone
        self.stack_SED(sed_blocks)

    def stack_SED(self, sed_blocks):
        # Order by incerasing metallicity order
        metals = [sb.Z for sb in sed_blocks]
        isrt_metal = np.argsort(metals)
        self.SED = np.stack([sb.SEDs for sb in np.array(sed_blocks)[isrt_metal]], axis=0)
        self.metal_points = np.array(metals)[isrt_metal]


class Sed_block():
    def __init__(self, fn="", load=True):
        self.fn = fn
        self.n_ages = 0
        self.n_wavs = 0
        self.age_points = None
        self.wavelengths = None
        self.isochrone = None
        self.X=0
        self.Y=0
        self.Z=0
        self.SEDs = None
        self.arr1 = None
        self.arr2 = None
        if load:
            self.load_bc()

    def load_bc(self, fn=None):
        """
        Loads SED template of BC2003 data.
        Note that the data is given in a very inconvenient form. :(
        """
        if fn is None:
            fn = self.fn
        #if fn[-2:] == "gz":
        import gzip
        f = gzip.open(fn)
        #else:
        #    f = open(fn, "r")

        f.seek(0)
        l = f.readline()
        lsp = l.split()
        self.n_ages = int(lsp[0])
        ages_l = l.split()[1:]

        while len(ages_l) < self.n_ages:
            l = f.readline()
            ages_l.extend(l.split())

        self.age_points = np.array([float(aa) for aa in ages_l])

        # Unknown lines
        ll = f.readline().decode("utf-8")
        #self.int1, self. int2, self.int3 = [int(ii) for ii in ll.split()]
        f.readline()
        #f.readline()
        cnt = 0
        while cnt < 3:
            desc = f.readline().decode("utf-8")
            if "(" in desc:
                cnt += 1

        self.isochrone = desc.split(").")[0]+")"
        xyz = desc.split(").")[1]
        self.X = float(xyz.split(",")[0].split("=")[1])
        self.Y = float(xyz.split(",")[1].split("=")[1])
        self.Z = float(xyz.split(",")[2].split("=")[1])


        lsp = f.readline().split()
        self.n_wavs = int(lsp[0])
        wav_l = lsp[1:]
        while len(wav_l) < self.n_wavs:
            wav_l.extend(f.readline().split())

        self.wavelengths = np.array([float(aa) for aa in wav_l])

        all_l = []
        for ll in f.readlines():
            all_l.extend(ll.split())

        arr1 = []
        SED = []
        i = 0
        i_age = 0
        while i_age < self.n_ages:
            n_wav = int(all_l[i])
            i +=1
            SED.append(np.array([float(aa) for aa in all_l[i:i+n_wav]]))
            i += n_wav
            # Something else
            n_something = int(all_l[i])
            i +=1
            arr1.append(np.array([float(aa) for aa in all_l[i:i+n_something]]))
            i += n_something
            i_age += 1

        self.SEDs = np.vstack(SED)
        self.arr1 = np.vstack(arr1)

        #i = 1536834
        #arr2 = np.zeros((10, n_ages))
        arr2=[]
        while i < len(all_l):
            na = all_l[i]
            i +=1
            arr2.append(np.array([float(aa) for aa in all_l[i:i+self.n_ages]]))
            i += self.n_ages

        self.arr2 = np.vstack(arr2)

        f.close()
        

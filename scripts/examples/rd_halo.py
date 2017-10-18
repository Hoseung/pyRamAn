
if __name__=="__main__":
    import tree.halomodule as hmo
    import matplotlib.pyplot as plt
    import numpy as np

    hcat = hmo.Halo(base="./29172/", nout=187, is_gal=False)
    plt.hist(hcat.data["mvir"])
    plt.suptitle("Halo mass")
    plt.show()


    # Tree
    
    



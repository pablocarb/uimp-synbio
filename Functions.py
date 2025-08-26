import numpy as np

def limonene_production_dynamics(y, t, rbs1_strength, rbs2_strength):
    """
    Simulates the limonene production model based on RBS strengths.
    
    Parameters:
    rbs1_strength (float): Strength of the first RBS.
    rbs2_strength (float): Strength of the second RBS.
    time (array): Time points for simulation.
    
    Returns:
    array: Limonene production values over time.
    """
    # Placeholder for the actual model implementation
    # This should include the ODEs and numerical integration logic

    gppsm, limsm, gpps, lims, gpp, lim = y  # Unpack the state variables

    ipp = 10  # Initial concentration of intermediate product, can be adjusted
    
    # Define parameters for the model
    Km1 = 0.1  # Transcription rate for first step
    Kdm1 = 0.01  # Degradation rate for first step
    Km2 = 0.1  # Transcription rate for second step
    Kdm2 = 0.01  # Degradation rate for second step
    Kr1 = np.arctan(rbs1_strength/20)  # RBS strength for first gene
    Kr1 = np.tanh(rbs1_strength/50)  # RBS strength for first gene
    Kd1 = 0.01  # Degradation rate for first gene product
    Kr2 = np.arctan(rbs2_strength/20)  # RBS strength for second gene
    Kr2 = np.tanh(rbs2_strength/50)  # RBS strength for first gene
    Kd2 = 0.01  # Degradation rate for second gene product
    
    Kcat = [0.01, 0.5]  # Catalytic constants for the reactions
    Km = [100, 1]  # Michaelis constants for the reactions
    Kd_gpp = 0.01  # Degradation rate for GPP
    Kd_lim = 0.01  # Degradation rate for limonene
    Kn1 = 1
    Kn2 = 1
    
    # Assuming a simple linear relationship for demonstration purposes  
    # In practice, this should be replaced with the actual ODE solution
    dgppsm = Km1 - Kdm1 * gppsm 
    dlimsm = Km2 - Kdm2 * limsm 
    dgpps = Kr1 * gppsm/(gppsm+Kn1) - Kd1 * gpps
    dlims = Kr2 * limsm/(limsm+Kn2) - Kd2 * lims  

    dgpp = Kcat[0] * gpps * ipp/(Km[0]+ipp) -  Kd_gpp * gpp
    dlim = Kcat[1] * lims * gpp/(Km[1]+gpp) -  Kd_lim * lim
    
    # Example: Limonene production is proportional to the product of RBS strengths
    # and time, simulating a simple growth model
    return [dgppsm, dlimsm, dgpps, dlims, dgpp, dlim]


def get_kmers(sequence, k=2):
        """
        Generate k-mers from a given sequence.
        
        Parameters:
        sequence (str): The input sequence.
        k (int): Length of the k-mers to generate.
        
        Returns:
        list: List of k-mers.
        """
        km = {}
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if kmer not in km:
                km[kmer] = []
            km[kmer].append(i)
        return km



def get_rbs_from_sequence(combinations, combinations1, combinations2):


    def get_kmers(sequence, k=2):
        """
        Generate k-mers from a given sequence.
        
        Parameters:
        sequence (str): The input sequence.
        k (int): Length of the k-mers to generate.
        
        Returns:
        list: List of k-mers.
        """
        km = {}
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if kmer not in km:
                km[kmer] = []
            km[kmer].append(i)
        return km
    
    k1 = set()
    k2 = set()
    k3 = set()

    for i in range(len(combinations)):
        seq = ''.join(combinations[i])
        kmers1 = get_kmers(seq, k=1)
        k1 |= set(kmers1)
        kmers2 = get_kmers(seq, k=2)
        k2 |= set(kmers2)
        kmers3 = get_kmers(seq, k=3)
        k3 |= set(kmers3)

    rbs1 = {}
    for k in k1:
        rbs1[k] = np.random.uniform()
    rbs2 = {}
    for k in k2:
        rbs2[k] = np.random.uniform()
    rbs3 = {}
    for k in k3:
        rbs3[k] = np.random.uniform()

    rbss1 = {}
    for i in range(len(combinations1)):
        seq = ''.join(combinations1[i])
        rbss1[combinations1[i]] = 1
        kmers1 = get_kmers(seq, k=1)
        for j in kmers1:
            for pos1 in kmers1[j]:
                rbss1[combinations1[i]] += 10*rbs1[j]
        kmers2 = get_kmers(seq, k=2)
        for j in kmers2:
            rbss1[combinations1[i]] += 0*rbs2[j]
        kmers3 = get_kmers(seq, k=3)
        for j in kmers3:
            rbss1[combinations1[i]] += 0*rbs3[j]

    rbss2 = {}
    for i in range(len(combinations2)):
        seq = ''.join(combinations2[i])
        rbss2[combinations2[i]] = 1
        kmers1 = get_kmers(seq, k=1)
        for j in kmers1:
            for pos2 in kmers1[j]:
                rbss2[combinations2[i]] += 10*rbs1[j] #*pos2
        kmers2 = get_kmers(seq, k=2)
        for j in kmers2:
            for pos2 in kmers2[j]:
                rbss2[combinations2[i]] += 0*rbs2[j] #*pos2**2
        kmers3 = get_kmers(seq, k=3)
        for j in kmers3:
            for pos2 in kmers3[j]:
                rbss2[combinations2[i]] += 0*rbs3[j]
    return rbss1, rbss2


# Create dummy input for the training set
def one_hot_encoding(mylib):
    Xd = []
    vn = ['A','C','G']
    for construct in mylib:
        v = []
        for pos in construct:
            for n in pos:         
                for i in np.arange(len(vn)): 
                    if n == vn[i]:
                        v.append(1)
                    else:
                        v.append(0)
        Xd.append(v)
    return np.array(Xd)
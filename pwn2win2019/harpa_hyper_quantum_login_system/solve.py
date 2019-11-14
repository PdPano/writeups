#Too lazy for markdown

import numpy as np
import qiskit
from scipy.linalg import logm
from collections import Counter

# Main steps to solution
# Identify the problem uses the HHL algorithm
# Identify that many operators can be reduced to identity
# Use a stack to simplify the code
# Extract the A matrix from the relevant portion of the code
# Recover key using matrix A and OTP

# Not super reliable, but works for this case
stack = ['//init/n']
for l in open('chall.qasm').readlines():
    l = l.replace('q2','b').replace('q3','c') # Regs identified in rev
    if l == stack[-1]:
        stack.pop()
        continue
    if l.startswith('cu3(0.412312561231198') and stack[-1].startswith('cu3(0.412312561231198') and stack[-1][-13:] == l[-13:]:
        stack.pop()
        continue
    if l.startswith('cu1(1.315689436019364') and stack[-1].startswith('cu1(-1.315689436019364') and stack[-1][-13:] == l[-13:]:
        stack.pop()
        continue
    if l.startswith('cu1(-pi)') and stack[-1].startswith('cu1(pi)') and stack[-1][-13:] == l[-13:]:
        stack.pop()
        continue
    stack.append(l)

def get_circuit(code):
    if qiskit.__version__=='0.6.1': # My dev env has a different version...
        c = qiskit.load_qasm_string(code)
    else:
        c = qiskit.QuantumCircuit.from_qasm_str(code)
    return c

def get_unitary(code):
    if isinstance(code,str):
        code = get_circuit(code)
    un = qiskit.Aer.get_backend('unitary_simulator')
    return qiskit.execute(code,un).result().get_unitary()

def get_A(U):
    return (-1.j*logm(U[1::2,1::2])).real/(2*np.pi)*2**8


# Must use the first iteration to extract matrix A.
# There is a condition that |max(eigv(A))| < 2**8 for the circuit to be valid,
# and in the first iteration we get A/2**8 in the exponent.
# If we don't use the first iteration, the matrix may have eigv > 1 and be impossible
# to invert.
def extract_code(code):
    return code[109:226] #can use code in https://github.com/nelimee/quantum-hhl-4x4/tree/master/hhl4x4 as a reference to figure the correct boundary

encoding_dict = { c:'{:06b}'.format(v) for v,c in enumerate('ABCDEFGHIJKLMNOPQRSTUVXWYZabcdefghijklmnopqrstuvxwyz01234567{}-_')}
decoding_dict = { v:k for k,v in encoding_dict.items() }

def encode_string(s):
    return ''.join(encoding_dict[c] for c in s )

def decode_string(s):
    dec = ''
    for i in range(0,len(s),6):
        dec += decoding_dict.get( s[i:i+6], '')
    return dec


new_code = '''OPENQASM 2.0;
              include "qelib1.inc";
              qreg c[1];
              qreg b[6];
'''           +''.join(extract_code(stack))
print(new_code)

print(np.linalg.norm(np.eye(64) - np.round(get_unitary(new_code)[::2,::2],3)))
print(np.round(get_unitary(new_code)[1::2,1::2],3).real)
A = get_A(get_unitary(new_code))
A = A[:16,:16] #Only need first 16x16 block

def get_otp():
    otp0 = list( list(v) for v in np.loadtxt('./otp/otp0.csv',delimiter=','))
    otp1 = list( list(v) for v in np.loadtxt('./otp/otp1.csv',delimiter=','))
    return [otp0,otp1]
otp = get_otp()
otp = np.array(otp)


# In a linear system we usually have A and b and want to solve for x
# A x = b
# 
# In this problem, we are given A and x, and must determine which
# b is to be used as an input, such that we get x as an answer.
#
# A quantum state of N qubits is encoded as a vector of 2**N entries
# So, for a 2-qubit state we would have:
# |psi> = a*|00> + b*|01> + c*|10> + d*|11> => (a,b,c,d)
#
# We have to check int(target[0][-NBITS-2:-1], 2) == int(kn)
# So, if we have kn=0, we want our state to maximize the chances
# of measuring |00>, and our target state should be (1,0,0,0)
#
# It is theoreticaly possible to use other states as a target and still
# get a match, but to get it right ~30 times per x and on all 16 columns
# it is not feasable in our lifetimes
#
# Now, given x corresponding to kn, we can determine
# that b is the kn'th column of the matrix A


# otp is related to the actual state up to a constant factor.
# Compute ratios to determine correct scaling factor
all_ratios = np.round(list((A/(otp[0]+1e-30)).ravel())+list((A/(otp[1]+1e-30)).ravel()),decimals=6)
counter = Counter(all_ratios)
print(f'Most frequent ratios: {counter.most_common(3)}')
otp_mul = otp * counter.most_common(2)[1][0]

#up to three decimal places is enough. Could use something like |a-b|<eps instead
A_approx = np.round(A,decimals=3)
otp_approx = np.round(otp_mul,decimals=3)

#inspect matrix and otp to determine the correct value here
#the extracted matrix A could be off by (128*j*np.eye(16))
np.fill_diagonal(A_approx,-32)


rec_keys = np.zeros_like(A_approx,dtype=np.int)
rec_keys[np.where(A_approx == otp_approx[1])] = 1

# Decoding what we got so far
print(decode_string(''.join([str(v) for v in rec_keys.T.ravel()])))

# It almost looks like a flag, but reversed
print(decode_string(''.join([str(v) for v in rec_keys.T.ravel()]))[::-1])

# Almost there, some bits are still flipped.
# otp is undecidable in 4 places
undecidable_inds = np.where(otp_approx[0]==otp_approx[1])
unpack_inds = list(zip(*undecidable_inds))
num_undec = len(unpack_inds)

print(unpack_inds)
# Brute force last 4 undecidable bits
for i in range(2**num_undec):
    c_rec_keys = rec_keys.copy()
    for j in range(num_undec):
        if i & (2**j):
            c_rec_keys[unpack_inds[j]] = 0
        orig_enc_string = ''.join([str(v) for v in c_rec_keys.T.ravel()])
    print(decode_string(orig_enc_string)[::-1])

# At last:
# CTF-BR{N0_mor3_LIn3aR_SysTeMs_to_SOlvePlz}



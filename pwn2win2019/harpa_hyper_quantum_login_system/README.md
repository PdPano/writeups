# Pwn2Win CTF 2019: HARPA Hyper quantum Login system
**DISCLAIMER:** *I'm the chall author, I just wrote the solution as if I was playing the CTF!*
**Category:** Rev
**Points:** 500
**Solves:** 0

> HARPA implemented quantum authentication systems for their machines. We found this prototype and believe it is possible to recover the password. We need your help.

## Hints
> Even quantum Systems are just Linear operators in disguise.
> Quantum Information algorithms can be implemented following any convention for the order of qubits in the registers.

## Write-up

### Rev
We are given 4 files: `control.py`, `chall.qasm`, `otp/otp{0,1}.csv`. After
reversing the `control.py` we conclude that we have partial control of 
the input to the quantum circuit in `chall.qasm`. Our input string is decoded
in a Base64 variant and the bits are used to choose elements from the matrices
stored in `otp/otp{0,1}.csv`.

```python
def encode_string(s):
    return ''.join(encoding_dict[c] for c in s )

def gen_keys(enc_s):
    if len(enc_s)<BLOCK_EL: enc_s +='0'*(BLOCK_EL-len(enc_s))
    return [enc_s[i:i+BLOCK_SIZE] for i in range(0,BLOCK_EL,BLOCK_SIZE)]

def get_states(keys,otp):
    return [[otp[int(i)][j][kn] for j,i in enumerate(k)] for kn,k in enumerate(keys)]
```

These selected elements were broken in 16 pieces, and now
each piece is appended with an unknown padding and the results are passed as input
to the circuit via the generated file `converted_input.inc`.

```python
def initial_state(s):
    res = np.array(s)
    return res/np.linalg.norm(res)

def setup_state(vec):
    out = ( 'initialize('
            + ','.join('{:.9f}'.format(v) for v in vec)
            + ')'
            + ','.join('b[{}]'.format(j) for j in range(NBITS))
            + ';\n'
            )
    return out

# main

keys = gen_keys(encode_string(pwd))
states = get_states(gen_keys(encode_string(pwd)), otp)

for kn,padding,state in zip(range(len(states)),my_secrets.padding, states):
    with open("converted_input.inc","w") as f:
        print(setup_state(initial_state(state + padding)))
        f.write(setup_state(initial_state(state + padding)))
```

The result of each run is filtered on one of the output values, it is checked that there are
counts for only one state, the resulting state is converted to an integer and compared
to the loop iteration counter. If we pass all tests in all rounds, we can login

```python
# In the main loop
target = []
while len(target)==0:
    results = qiskit.execute(circ, backend, shots=30).results()
    counts = results.get_counts(circ)
    target = list(filter(lambda s: s[-1]=='1', counts.keys()))

if len(target)==1 and int(target[0][-NBITS-2:-1], 2) == int(kn):
    count+=1

# Final check
if count==BLOCK_SIZE:
    print("Login successful!")
```

### But, what is it actually doing?
There are many examples of quantum algorithms, but not that many have been fully implemented
in qasm. Also, the capitalization in the title spells HHL... After some research, we
determine that the circuit implements the HHL algorithm for solving linear systems of the form
`A*x=b`. We are finally able to understand that the circuit is using the state generated
from our password and using it as the `b`. The matrix `A` must be extracted from the circuit,
`x` is the `target` and we must determine which `b` will make the test
`int(target[0][-NBITS-2:-1], 2) == int(kn)` pass.

We must also satisfy `len(target)==1`, which means our state must have only one component.
To understand why, we need some quantum mechanics.
A quantum state of N qubits is encoded as a vector of 2^N entries.
So, for a 2-qubit state we would have:
`|psi> = a*|00> + b*|01> + c*|10> + d*|11> => (a,b,c,d)`.
The sum of the modulus-squared of all coefficients of a state `|psi>` must
equal one, and the probability of measuring a certain state is the
modulus square of it's coefficient.
So, if we have `kn=0`, we want our state to maximize the chances
of measuring `|00>`, our target state should be (1,0,0,0).
Finally, if we multiply `A` by a vector with only one non-zero entry
we get a column of A as a result, and that is our `b` for that round!


### Extracting A
Looking at the qasm code, some parts look random, and some parts seem to follow a pattern.
It was possible to identify that the random parts can be reduced to no-ops, leaving a much
more treatable circuit. Also, we could find an
[implementation of the algorithm](https://github.com/nelimee/quantum-hhl-4x4/tree/master/hhl4x4) 
that follows the same patterns present in the processed circuit.

Using the found implementation code as a reference,
we extracted the matrix `A` following the procedure suggested there.
```python
def get_unitary(code):
    if isinstance(code,str):
        code = get_circuit(code)
    un = qiskit.Aer.get_backend('unitary_simulator')
    return qiskit.execute(code,un).result().get_unitary()

def get_A(U):
    return (-1.j*logm(U[1::2,1::2])).real/(2*np.pi)*2**8

A = get_A(get_unitary('relevant code'))
A = A[:16,:16] #Only need first 16x16 block
```

### Recovering the password
The state we construct from our password is equal to the columns in matrix A
up to a multiplicative constant. We can find the most frequent ratios
to recover the right constant and retrieve the bit pattern. (Note: in principle,
it is possible to have different constants for each column. For reasons I won't go into,
I decided to use a single constant when building the chall).

```python
all_ratios = np.round(list((A/(otp[0]+1e-30)).ravel())+list((A/(otp[1]+1e-30)).ravel()),decimals=6)
counter = Counter(all_ratios)
otp_mul = otp * counter.most_common(2)[1][0] #Most frequent is zero
A_approx = np.round(A,decimals=3)
otp_approx = np.round(otp_mul,decimals=3)
np.fill_diagonal(A_approx,-32) # Euler's identity can mess up the diagonal. Use otp values to determine the correct one

rec_keys = np.zeros_like(A_approx,dtype=np.int)
rec_keys[np.where(A_approx == otp_approx[1])] = 1
decoding_dict = { v:k for k,v in encoding_dict.items() }

def decode_string(s):
    dec = ''
    for i in range(0,len(s),6):
        dec += decoding_dict.get( s[i:i+6], '')
    return dec

print(decode_string(''.join([str(v) for v in rec_keys.T.ravel()])))
```

This prints `}zlPevlOS_ot_sMeTsyS_Ra3nML__rom_0t{RB-HTC`, which when reversed
gives `CTH-BR{t0_mor__LMn3aR_SysTeMs_to_SOlvePlz}`. This is very close to a flag
and we could actually brute-force all the letters that seem to be wrong, but there is
a better way. Turns out the OTP is undecidable in four locations (`otp[0][i,j]==otp[1][i,j]` for 4 pairs `i,j`)
and we only need to brute-force these 4 bits. There is only one string that looks like a flag: `CTF-BR{N0_mor3_LIn3aR_SysTeMs_to_SOlvePlz}`


### Additional comment

Right after the CTF, I talked with the teams that were trying to solve the chall to see how close they were.
They had the idea of removing the control bit to process the unitary matrix. While developing the chall I though this could be done,
but never tried to implement it as I already had another solver.

I was really bothered after seeing that this solution was not working and I tried to find something wrong with it.
At first, I thought that the operator `tdg c[0];` could not be ignored as they had done, and that's what I told them, but I was wrong.
In fact, they transformed `ccx c[0],q2[0],q2[1];` in `cx q2[0],q2[1];` and this was the real issue. Originally, if `c[0]` is zero
`q2[1]` will not flip irrespective of `q2[0]`'s state, but in the modified version this is no longer the case!

In the end, I kept learning things even after the CTF!

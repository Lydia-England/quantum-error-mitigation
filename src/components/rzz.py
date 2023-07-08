#--- Build and Return ZZ Rotation Gate ---#
def get_rzz():
    qr_zz = QuantumRegister(2, 'qc_zz')
    qc_zz = QuantumCircuit(qr_zz)
    qc_zz.append(SDG,[0])
    qc_zz.append(SDG,[1])
    qc_zz.append(SY,[1])
    qc_zz.cx(qr_zz[0],qr_zz[1])
    qc_zz.append(SYDG,[1])
    # add a global phase to correct for global phase introduced
    qc_zz.global_phase += pi/4
    return qc_rzz


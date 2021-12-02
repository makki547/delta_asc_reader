import numpy as np
import nmrglue as ng
import re

def _read_hdr(hdr_filename):

    string_pattern = re.compile(r'"(.*)"')
    integer_pattern = re.compile(r'')

    val_patterns = [
        (re.compile(r'"(TRUE|FALSE)"'), bool, False),
        (re.compile(r'"(.*)"'), str, False),
        (re.compile(r'(\d+\.?\d*)\[(.+)\]'), float, True),
        (re.compile(r'(\d+)'), int, False),
        (re.compile(r'(\d+\.\d*)'), float, False)
    ]

    
    header = {}
    with open(hdr_filename) as fin:

        for line in fin:

            elems = line.strip().split()

            if len(elems) < 2:
                continue


            val = ''.join(elems[1:])

            for pattern, conv, has_unit in val_patterns:


                match = pattern.match(val)

                if match:
                    val = conv(match.group(1))

                    if has_unit:
                        unit = match.group(2)
                    else:
                        unit = None

                    header[elems[0]] = (val, unit)
                    break


    return header


def _process_si_prefix(val, unit_with_si):

    si_factor = {'p':1.0E-12, 'n': 1.0E-9, 'u': 1.0E-6, 'm':1.0E-3,
                 'k': 1.0E+3, 'M': 1.0E+6, 'G':1.0E+12, 'T':1.0E+15}
    si_unit_pattern = re.compile(r'(p|n|u|m|k|M|G|T)(.+)')

    if unit_with_si is None:
        return val, unit_with_si
    
    #"ppm" matches the pattern, then escape 
    if unit_with_si == 'ppm':
        return val, 'ppm'

    match = si_unit_pattern.match(unit_with_si)
    if match:
        return val * si_factor[match.group(1)], match.group(2)

    else:
        return val, unit_with_si


    
def _read_complex1d(asc_file_name, size):

    data = np.zeros(size, dtype = np.complex64)

    with open(asc_file_name) as fin:

        idx = 0
        for i, line in enumerate(fin):


            if i == 0: continue

            elems = line.strip().split()
            if len(elems) != 3:
                raise
            
            data[i - 1] = float(elems[1]) + 1.0j * float(elems[2])
                

    return data

def _read_real1d(asc_file_name, size):

    data = np.zeros(size, dtype = np.float32)

    with open(asc_file_name) as fin:

        idx = 0
        for i, line in enumerate(fin):


            if i == 0: continue

            elems = line.strip().split()
            if len(elems) != 2:
                raise
            
            data[i - 1] = float(elems[1])
                

    return data

def _read_complex2d(asc_file_name, size):

    data = np.zeros(size, dtype = np.complex64)

    with open(asc_file_name) as fin:

        x_idx = 0
        y_idx = 0
        for i, line in enumerate(fin):


            if i == 0: continue

            elems = line.strip().split()
            if len(elems) != 4:
                raise
            
            data[y_idx, x_idx] = float(elems[2]) + 1.0j * float(elems[3])

            x_idx += 1
            if x_idx == size[1]:
                x_idx = 0
                y_idx += 1
                

    return data

def _read_real2d(asc_file_name, size):

    data = np.zeros(size, dtype = np.float32)

    with open(asc_file_name) as fin:

        x_idx = 0
        y_idx = 0
        for i, line in enumerate(fin):

            if i == 0: continue

            elems = line.strip().split()
            if len(elems) != 3:
                raise
            
            data[y_idx, x_idx] = float(elems[2])

            x_idx += 1
            if x_idx == size[1]:
                x_idx = 0
                y_idx += 1            

                

    return data

def _interleave_2d(data):

    il_data = np.zeros(data.shape, dtype = data.dtype)
    il_data[0::2,:] = np.conjugate(data[:data.shape[0]//2,:])
    il_data[1::2,:] = -1.0 * np.conjugate(data[data.shape[0]//2:,:])
    return il_data

def _delete_2d_sine_modulation(data):

    cos_data = np.zeros((data.shape[0] // 2, data.shape[1]), dtype = data.dtype)
    cos_data = data[0:data.shape[0]//2,:]

    return cos_data
            


# Read .asc and its related .hdr files made by Delta
# and output dic and data of nmrglue.pipe
def read(asc_filename):

    filename_pattern = re.compile(r'(.+)\.asc')

    match = filename_pattern.match(asc_filename)
    if match is None:
        raise


    header = _read_hdr(match.group(1) + '.hdr')


    ndim = header['dimensions'][0]
    if ndim > 3:
        #raise here exception
        return None, None

    udic = {'ndim': ndim}
    udic.update({i:{} for i in range(ndim)})
    dim_idx = {'x': ndim - 1, 'y': ndim - 2}
    
    x_freq, _ = _process_si_prefix(*header['x_freq'])

    x_idx = dim_idx['x']
    udic[x_idx]['obs'] = x_freq / 1.0E+6
    udic[x_idx]['car'] = x_freq / 1.0E+6 * header['x_offset'][0]
    udic[x_idx]['sw'] = _process_si_prefix(*header['x_sweep'])[0]
    udic[x_idx]['label'] = header['x_domain'][0]
    udic[x_idx]['complex'] = (header['x_format'][0] == 'COMPLEX')
    udic[x_idx]['freq'] = (header['x_start'][1] == 'ppm')
    udic[x_idx]['time'] = not udic[x_idx]['freq']
    udic[x_idx]['size'] = header['x_curr_points'][0]
    udic[x_idx]['encoding'] = 'complex'
        

    if ndim == 2:
        y_idx = dim_idx['y']
        y_freq, _ = _process_si_prefix(*header['y_freq'])
        udic[y_idx]['obs'] = y_freq / 1.0E+6
        udic[y_idx]['car'] = y_freq / 1.0E+6 * header['y_offset'][0]
        udic[y_idx]['sw'] = _process_si_prefix(*header['y_sweep'])[0]
        udic[y_idx]['label'] = header.get('y_domain', ('None',))[0]
        udic[y_idx]['complex'] = (header['y_format'][0] == 'COMPLEX')
        udic[y_idx]['freq'] = (header['y_start'][1] == 'ppm')
        udic[y_idx]['time'] = not udic[y_idx]['freq']
        udic[y_idx]['size'] = header['y_curr_points'][0]
        if udic[y_idx]['complex']:
            udic[y_idx]['size'] *= 2
            
        udic[y_idx]['encoding'] = 'states'

    if ndim == 1:

        if udic[x_idx]['complex']:
            data = _read_complex1d(asc_filename, udic[x_idx]['size'])
        else:
            data = _read_real1d(asc_filename, udic[x_idx]['size'])
            
    elif ndim == 2:

        if udic[x_idx]['complex']:
            data = _read_complex2d(asc_filename, (udic[y_idx]['size'], udic[x_idx]['size']) )
        else:
            data = _read_real2d(asc_filename, (udic[y_idx]['size'], udic[x_idx]['size']) )


        if udic[y_idx]['complex'] and udic[y_idx]['time']:
            data = _interleave_2d(data)
        if udic[y_idx]['complex'] and udic[y_idx]['freq']:
            udic[y_idx]['size'] //= 2
            udic[y_idx]['complex'] = False
            data = _delete_2d_sine_modulation(data)
            
        if udic[y_idx]['label'] == 'None':
            udic[y_idx]['obs'] = 1.0


    pipe_dic = ng.pipe.create_dic(udic)

    return pipe_dic, data


def get_array_list(asc_filename, dim = 'y'):

    filename_pattern = re.compile(r'(.+)\.asc')

    match = filename_pattern.match(asc_filename)
    if match is None:
        raise


    header = _read_hdr(match.group(1) + '.hdr')


if __name__ == '__main__':
    import sys

    #read asc file
    dic, data = read(sys.argv[1])
    

    ng.pipe.write('test2.ft2', dic, data, overwrite= True)

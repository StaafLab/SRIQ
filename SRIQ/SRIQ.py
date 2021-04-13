import os

def runSRIQ(data, cutOff = None, permutations = 100, iterations = 3, minBagSize=900, minClusterSize = 0):
    #data: string with the directory to the csv file which you want to run
    #cutOff: list of the cutoffs
    output = ''
    resources = 'software/VRLA/resources/test.properties'
    if cutOff is None: cutOff = [0.48, 0.50]
    with open(resources) as file:
        for line in file.readlines():
            if 'inFileName' in line:
                output += f'inFileName={data}\n'
            elif 'distCutOff' in line:
                if isinstance(cutOff, list): 
                    output += 'distCutOff={}\n'.format(", ".join([str(x) for x in cutOff]))
                else:
                    output += f'distCutOff={str(cutOff)}\n'
            elif 'permutations' in line:
                output += f'permutations={permutations}\n'
            elif 'minClusterSize' in line:
                output += f'minClusterSize={minClusterSize}\n'
            elif 'minBagSize' in line:
                output += f'minBagSize={minBagSize}\n'
            elif 'iterations' in line:
                output += f'iterations={iterations}\n'
            else:
                output += line
    with open(resources, 'w') as file:
        file.write(output)
    print('running...')
    bashCommand = "cd ../software/VRLA && java -jar -Xmx4g VRLA.jar"
    os.system(bashCommand)
    print('Done!')
    
runSRIQ(data = 'test_mc_log2var(80)')
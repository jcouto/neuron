function [neurons,input,info] = loadCommonInputFile(filename,verbose)
% [neurons,input,info] = loadCommonInputFile(filename,verbose=0)
if ~exist('verbose','var'),verbose=0;end
input = struct([]);

try
    input.spikes = hdf5read(filename, '/Input/PresynapticSpikes');
catch
    printError('/Input/PresynapticSpikes',verbose)
end

try
    props = hdf5read(filename, '/Input/Properties');
    for k=1:length(props.MemberNames)
        input.(props.MemberNames{k}) = props.Data{k};
    end
catch
    printError('/Input/Properties',verbose)
end


info = struct([]);
try
    props = hdf5read(filename, '/Properties/Simulation');
    for k=1:length(props.MemberNames)
        info(1).(props.MemberNames{k}) = props.Data{k};
    end
catch
    printError('/Input/Simulation',verbose)
end

%neurons = repmat(struct([]), [info.nNeurons,1]);
nrninfo = h5info(filename,'/Neurons');
N = length(nrninfo.Groups);
neurons = repmat(struct([]), [N,1]);
for k=1:N
    neurons(k).spikes = hdf5read(filename, sprintf('/Neurons/ID_%d/Spikes', k-1));
    try
        neurons(k).pert = hdf5read(filename, sprintf('/Neurons/ID_%d/Perturbations', k-1));
    catch
        printError(sprintf('/Neurons/ID_%d/Perturbations', k-1),verbose)
    end
end
try
    neurons(k).v = hdf5read(filename, sprintf('/Neurons/ID_%d/Voltage', k-1));
catch
    printError(sprintf('/Neurons/ID_%d/Voltage', k-1),verbose)
end
end

function printError(x,verbose)
    if(verbose)
        fprintf(1,'loadCommonInputFile: %s not found.\n',x)
    end
end
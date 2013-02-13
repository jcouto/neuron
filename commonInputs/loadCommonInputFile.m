function [neurons,input,info] = loadCommonInputFile(filename)
% [neurons,input,info] = loadCommonInputFile(filename)

input.spikes = hdf5read(filename, '/Input/PresynapticSpikes');
props = hdf5read(filename, '/Input/Properties');
for k=1:length(props.MemberNames)
    input.(props.MemberNames{k}) = props.Data{k};
end

info = struct([]);
props = hdf5read(filename, '/Properties/Simulation');
for k=1:length(props.MemberNames)
    info(1).(props.MemberNames{k}) = props.Data{k};
end

neurons = repmat(struct([]), [info.nNeurons,1]);
for k=1:info.nNeurons
    neurons(k).spikes = hdf5read(filename, sprintf('/Neurons/ID_%d/Spikes', k-1));
    try
        neurons(k).pert = hdf5read(filename, sprintf('/Neurons/ID_%d/Perturbations', k-1));
    catch
    end
    try
        neurons(k).v = hdf5read(filename, sprintf('/Neurons/ID_%d/Voltage', k-1));
    catch
    end
end


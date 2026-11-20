function[] = vars2GPU()
    assignin('base','HxIncList',gpuArray(evalin('base','HxIncList')));
    assignin('base','EzIncList',gpuArray(evalin('base','EzIncList')));
    
    assignin('base','EzMatrix',gpuArray(evalin('base','EzMatrix')));
    assignin('base','KEzEz',gpuArray(evalin('base','KEzEz')));
    assignin('base','KEzHy',gpuArray(evalin('base','KEzHy')));
    assignin('base','KEzHx',gpuArray(evalin('base','KEzHx')));

    assignin('base','EzxPML',gpuArray(evalin('base','EzxPML')));
    assignin('base','PEzxEzx',gpuArray(evalin('base','PEzxEzx')));
    assignin('base','PEzxHy',gpuArray(evalin('base','PEzxHy')));
    assignin('base','EzyPML',gpuArray(evalin('base','EzyPML')));
    assignin('base','PEzyEzy',gpuArray(evalin('base','PEzyEzy')));
    assignin('base','PEzyHx',gpuArray(evalin('base','PEzyHx')));
    
    assignin('base','HxMatrix',gpuArray(evalin('base','HxMatrix')));
    assignin('base','KHxHx',gpuArray(evalin('base','KHxHx')));
    assignin('base','KHxEz',gpuArray(evalin('base','KHxEz')));

    assignin('base','HyMatrix',gpuArray(evalin('base','HyMatrix')));
    assignin('base','KHyHy',gpuArray(evalin('base','KHyHy')));
    assignin('base','KHyEz',gpuArray(evalin('base','KHyEz')));

    assignin('base','EzMonList',gpuArray(evalin('base','EzMonList')));
    assignin('base','HxMonList',gpuArray(evalin('base','HxMonList')));
    assignin('base','HyMonList',gpuArray(evalin('base','HyMonList')));
    assignin('base','totSCSList',gpuArray(evalin('base','totSCSList')));
end
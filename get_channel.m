function [prm, channel] = get_channel(prm)
% return a raytracing channel 
    
    rtchan = comm.RayTracingChannel(prm.rays{1}, prm.txSite, prm.rxSite);
    channel = rtchan;
end


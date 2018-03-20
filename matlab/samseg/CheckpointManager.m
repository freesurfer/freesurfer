classdef CheckpointManager < handle
   properties
      checkpoint_state
   end
   methods
      function obj = CheckpointManager()
         obj.checkpoint_state = {};
      end
      function save(obj, checkpoint_name, varnames)
        matlabDumpDir = getenv('MATLAB_DUMP_DIR');
        if ~isfield(obj.checkpoint_state, checkpoint_name)
            obj.checkpoint_state.(checkpoint_name) = 0;
        end
        obj.checkpoint_state.(checkpoint_name) = obj.checkpoint_state.(checkpoint_name) + 1;
        path =  [matlabDumpDir '/' checkpoint_name '_' num2str(obj.checkpoint_state.(checkpoint_name)) '.mat'];
        cmd = ['save ''' path ''' ' varnames];
        evalin('caller', cmd);
        disp(['Saving MATLAB dump to ' path]);
      end
   end
 end

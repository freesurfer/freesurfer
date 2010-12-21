-- A module to hold my Freesurfer related files

module( ..., package.seeall )

-- Returns the path to the scratch directory
function ScratchDir()
   return os.getenv( "SCRATCHDIR" )
end


-- Returns the location of the directory holding GCAM files
function GCAMdir()
   return ScratchDir().."/gcam/"
end

-- Returns the location of the directory holding MRI files
function MRIdir()
   return ScratchDir().."/mri/"
end

-- Returns the location of the registration directory
function Registrationdir()
   return ScratchDir().."/registration/"
end


-- Returns the location of the linear registration directory
function LinearRegistrationdir()
   return Registrationdir().."/linear/"
end

-- Returns the location of the mrilabel directory
function MRIlabelDir()
   return ScratchDir().."/mrilabel/"
end



-- Lists the files matching 'pattern' in the specified directory
function GetFiles( srcDir, pattern )
   print (srcDir)
   local lsPipe = io.popen( "ls "..srcDir )

   local fileList = {}

   local currFile
   for currFile in lsPipe:lines() do
      local trimName
      for trimName in currFile:gmatch( pattern ) do
	 table.insert( fileList, trimName )
      end
   end

   return fileList
end


-- Returns all the GCAM files
function AllGCAMfiles()
   return GetFiles( GCAMdir(), "(gcam%d+).nc" )
end


-- From lua-users wiki
-- Splits a string at the given pattern
function split(str, pat)
   local t = {}  -- NOTE: use {n = 0} in Lua-5.0
   local fpat = "(.-)" .. pat
   local last_end = 1
   local s, e, cap = str:find(fpat, 1)
   while s do
      if s ~= 1 or cap ~= "" then
	 table.insert(t,cap)
      end
      last_end = e+1
      s, e, cap = str:find(fpat, last_end)
   end
   if last_end <= #str then
      cap = str:sub(last_end)
      table.insert(t, cap)
   end
   return t
end
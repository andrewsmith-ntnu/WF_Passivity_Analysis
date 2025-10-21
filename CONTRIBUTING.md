Contributing

Thank you for contributing to this repository. Please follow these guidelines.

Branching
- Use feature branches: `feat/<short-description>`
- Use bugfix branches: `fix/<short-description>`

Commit messages
- Use imperative mood and be concise.
- Example: `Add docs generator and README`

Code style
- MATLAB functions: include a header comment block at the top with a short description, usage, inputs and outputs.

MATLAB function header template

% function_name  Short one-line description
%
%   [out1, out2] = function_name(in1, in2)
%
% Inputs
%   in1 - Description
%   in2 - Description
%
% Outputs
%   out1 - Description
%   out2 - Description
%
% Notes
%   Any notes or references.

Running the docs generator
- From MATLAB, run:
  cd(fullfile(fileparts(mfilename('fullpath'))));
  docs/generate_docs

Style notes
- Keep functions small and documented.
- Add unit tests when modifying core algorithms where possible.

Thank you!

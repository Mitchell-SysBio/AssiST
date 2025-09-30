function [reClassWell, reClass] = inputDialog2(myColorMap, numClasses)
% inputDialog2 — Positioned UI for manual reclassification with a color key
% Summary:
%   Shows a COLOR KEY (swatches + labels) for the 7 classes, then collects
%   well IDs and class IDs in a modal window. Live-validates inputs and
%   DISABLES the OK button when there is any error. Returns empty outputs
%   on cancel/close, and also when OK on both-empty inputs.
%
% Inputs:
%   promptTitle   char/string   Window title (optional; default: 'Reclassify wells').
%   position      1x4 double    Normalized [left bottom width height]
%                               (optional; default [.7 .5 .3 .4]).
%
% Outputs:
%   reClassWell   string row    Validated well IDs (e.g., ["A1","B3","H12"]) or [].
%   reClass       double row    Validated class IDs (1..7), aligned with wells, or [].

        promptTitle = 'Reclassify wells';
        position = [.7 .5 .3 .5]; % normalized

    % Default outputs (in case of cancel/close)
    reClassWell = string.empty(1,0);
    reClass     = [];

    % === Color/key data (matches rest of pipeline) =======================
    classNames = string(1:numClasses);
    % classNames = [ ...
    %     "1_bubbles"       ; ...
    %     "2_empty"         ; ...
    %     "3_pellet_tiny"   ; ...
    %     "4_cloud_small"   ; ...
    %     "5_pellet_small"  ; ...
    %     "6_cloud"         ; ...
    %     "7_pellet"        ];
    % % RGB rows in [0..1], consistent with inferPlatePhenotypes/plots
    % myColorMap = [ ...
    %     0.0000 0.0000 0.0000;  ... % black (bubbles)
    %     0.7882 0.0941 0.2902;  ... % red (empty)
    %     1.0000 0.4588 0.5608;  ... % pink (pellet_tiny)
    %     0.6549 0.7882 0.3412;  ... % light green (cloud_small)
    %     0.5647 0.8784 0.9373;  ... % light blue (pellet_small)
    %     0.0000 0.4471 0.0000;  ... % green (cloud)
    %     0.0000 0.4667 0.7137];     % blue  (pellet)

    % Build modal UI at requested position
    fig = uifigure('Name', promptTitle, 'Visible', 'off');
    fig.CloseRequestFcn = @onCancel;
    try
        fig.Units = 'normalized';
        fig.Position = position;
    catch
        % Fallback if normalized units unsupported (older versions)
    end

    % Top-level grid:
    % [ Color Key | Inputs/Validation/Buttons ]
    % Use 1 column; stack sections vertically for narrower dialogs
    gl = uigridlayout(fig, [3,1]);
    gl.RowHeight  = {'fit','fit','fit'};
    gl.Padding    = [12 12 12 12];
    gl.RowSpacing = 10;

    %% 1) COLOR KEY (swatches + labels)
    panelKey = uipanel(gl, 'Title', 'Color key (class → color)');
    keyGrid  = uigridlayout(panelKey, [7, 2]);
    keyGrid.ColumnWidth = {30, '1x'};
    keyGrid.RowHeight   = repmat({'fit'}, 1, 7);
    keyGrid.Padding     = [8 8 8 8];
    keyGrid.RowSpacing  = 4;
    keyGrid.ColumnSpacing=8;

    % Create 7 rows: swatch + label "k_name"
    for k = 1:numel(classNames)
        % Colored swatch (use a blank, bordered label as the swatch)
        sw = uilabel(keyGrid, 'Text', ' ');
        sw.BackgroundColor = myColorMap(k, :);
        %sw.BorderColor     = [0.3 0.3 0.3];
        sw.Layout.Row      = k;  sw.Layout.Column = 1;
        sw.Tooltip         = sprintf('RGB = [%.3f %.3f %.3f]', myColorMap(k,1), myColorMap(k,2), myColorMap(k,3));

        % Label text
        lbl = uilabel(keyGrid, 'Text', sprintf('%s', classNames(k)));
        lbl.Layout.Row      = k;  lbl.Layout.Column = 2;
    end

    %% 2) PROMPTS + FIELDS + INLINE MESSAGE
    inputPanel = uipanel(gl);
    form = uigridlayout(inputPanel, [5,1]);
    form.RowHeight  = {'fit','fit','fit','fit','fit'};
    form.RowSpacing = 6;
    form.Padding    = [0 0 0 0];

    % (Keep a short legend line for text-only clarity)
    uilabel(form, 'Text', 'Enter well IDs (A1–H12). Delimit with space/comma/semicolon:');
    wellsField = uieditfield(form, 'text');

    uilabel(form, 'Text', sprintf('Enter class IDs (1–%i). Same count as wells or a single value:', numClasses));
    classField = uieditfield(form, 'text');

    % Message area (inline errors/warnings)
    msg = uilabel(form, 'Text', '', 'FontColor', [0.85 0.33 0.10], 'WordWrap', 'on');

    %% 3) TIP + BUTTONS
    ctrlPanel = uipanel(gl);
    ctrl = uigridlayout(ctrlPanel, [1,3]);
    ctrl.ColumnWidth = {'1x','fit','fit'};

    tip = uilabel(ctrl, 'Text', 'Tip: leave both fields empty and press OK to skip reclassification.');
    tip.FontColor = [0.4 0.4 0.4];

    cancelBtn = uibutton(ctrl, 'Text', 'Cancel', 'ButtonPushedFcn', @onCancel);
    okBtn     = uibutton(ctrl, 'Text', 'OK',     'ButtonPushedFcn', @onOK, 'FontWeight', 'bold');

    % Live validation when user edits either field
    wellsField.ValueChangedFcn = @(~,~) validateForm();
    classField.ValueChangedFcn = @(~,~) validateForm();

    fig.Visible = 'on';
    validateForm();  % set initial OK state
    uiwait(fig);     % modal

    % If user pressed OK, data will be in fig.UserData
    if isvalid(fig) && ~isempty(fig.UserData) && isfield(fig.UserData,'status') && strcmp(fig.UserData.status,'ok')
        reClassWell = fig.UserData.reClassWell;
        reClass     = fig.UserData.reClass;
    end

    % Cleanup
    if isvalid(fig); delete(fig); end
    return

    %==== Nested callbacks and helpers ====================================

    function onCancel(~, ~)
        % Leave defaults (empties) and resume
        if isvalid(fig); fig.UserData = struct('status','cancel'); end
        uiresume(fig);
    end

    function onOK(~, ~)
        % Safety: run validation again. If invalid, do nothing (OK stays disabled anyway).
        [valid, wellsOut, classOut] = validateForm();
        if ~valid
            return
        end
        fig.UserData = struct('status','ok', 'reClassWell', wellsOut, 'reClass', classOut);
        uiresume(fig);
    end

    function [isValid, wellsOut, classOut] = validateForm()
        % Returns:
        %   isValid  — true if OK should be enabled
        %   wellsOut — string row (possibly empty)
        %   classOut — double row (possibly empty)

        % Normalize delimiters and case
        normDelims = @(s) regexprep(strtrim(s), '[,\t;]+', ' ');
        wellsStr = upper(normDelims(wellsField.Value));
        classStr = normDelims(classField.Value);

        % Tokenize
        wellsTokens = string(strsplit(wellsStr));
        wellsTokens(wellsTokens=="") = [];
        classTokens = strsplit(classStr);
        classTokens(strcmp(classTokens,"")) = [];

        % Allow both empty → valid, returns empties
        if isempty(wellsTokens) && isempty(classTokens)
            msg.Text = '';
            okBtn.Enable = 'on';
            isValid  = true;
            wellsOut = string.empty(1,0);
            classOut = [];
            return
        end

        % Validate wells
        if isempty(wellsTokens)
            msg.Text = 'Please enter well IDs (A1–H12) or leave both fields empty.';
            okBtn.Enable = 'off';
            isValid = false; wellsOut = []; classOut = [];
            return
        end
        isValidWell = @(tok) ~isempty(regexp(tok,'^[A-H](?:[1-9]|1[0-2])$','once'));
        validMask   = arrayfun(@(t) isValidWell(char(t)), wellsTokens);
        invalidWells = wellsTokens(~validMask);
        wellsTokens  = wellsTokens(validMask);

        if ~isempty(invalidWells)
            msg.Text = "Invalid well IDs: " + strjoin(cellstr(invalidWells), ', ');
            okBtn.Enable = 'off';
            isValid = false; wellsOut = []; classOut = [];
            return
        end

        % Validate classes
        if isempty(classTokens)
            msg.Text = sprintf('Please enter class IDs (1–%i) or leave both fields empty.',numClasses);
            okBtn.Enable = 'off';
            isValid = false; wellsOut = []; classOut = [];
            return
        end
        classVals = cellfun(@str2double, classTokens);
        validClassMask = ~isnan(classVals) & classVals>=1 & classVals<=7 & (mod(classVals,1)==0);
        invalidClassTokens = classTokens(~validClassMask);
        classVals = classVals(validClassMask);

        if ~isempty(invalidClassTokens)
            msg.Text = sprintf("Invalid class IDs (must be integers 1..%i): ",numClasses) + strjoin(invalidClassTokens, ', ');
            okBtn.Enable = 'off';
            isValid = false; wellsOut = []; classOut = [];
            return
        end
        if isempty(classVals)
            msg.Text = sprintf('No valid class IDs provided (must be integers 1..%i).',numClasses);
            okBtn.Enable = 'off';
            isValid = false; wellsOut = []; classOut = [];
            return
        end

        % Count logic: either single class (broadcast) OR counts must match
        nW = numel(wellsTokens);
        nC = numel(classVals);
        if ~(nC == 1 || nC == nW)
            msg.Text = sprintf('Count mismatch: wells=%d, classes=%d. Provide one class for all wells or the same count.', nW, nC);
            okBtn.Enable = 'off';
            isValid = false; wellsOut = []; classOut = [];
            return
        end

        % If we get here, inputs are valid → compute outputs
        if nC == 1 && nW > 1
            classVals = repmat(classVals, 1, nW);
        end
        msg.Text = '';            % clear any previous error
        okBtn.Enable = 'on';      % allow OK
        isValid  = true;
        wellsOut = wellsTokens(:).';
        classOut = double(classVals(:).');
    end
end

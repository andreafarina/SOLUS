current = cd;

TOAST_DIR = [];% path to TOAST
if ~isempty(TOAST_DIR)
    cd(TOAST_DIR)
    mtoast2_install
end

MATIMAGE_DIR  = [];%path to matImage
if ~isempty(TOAST_DIR)
    cd(MATIMAGE_DIR)
    installMatImage
end
cd(current)
clear current





from modelfunc import myerror
import yaml


def readModelParams(paramsFile, key=None):
    """Parse a model params dict from a yaml file. File can either have a
    single unnamed dict:
    x: abc
    y: xyz
    or several:
    key1:
        x: abc
        y: xyc
    key2:
        var1: 1
        var2: 5
    In the former case no key is need. In the latter case, if a specific dict
    is required, it should be specified using "key", otherwise a dict of dicts
    will be returned.
    Parameters
    ----------
    paramsFile : str
        yaml file with desired model params
    key: str [Optional]
        key to select which dict to return.

    Returns
    -------
    dict
        dict with model params with 'params' added as file it was read from.
    """
    if paramsFile is None:
        return {}
    try:
        with open(paramsFile) as fp:
            modelParams = yaml.load(fp, Loader=yaml.FullLoader)
            # print(modelParams)
            if key is not None:
                if key in modelParams:
                    modelParams = modelParams[key]
            # Force to use name of file actually being read
            modelParams['params'] = paramsFile
    except Exception:
        myerror(f'Could not open params file: {paramsFile}')
    return modelParams

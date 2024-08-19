import numpy as np
import torch
import torchvision.transforms.functional as F
import torchvision.transforms as T
from torchvision.models.optical_flow import raft_large
from torchvision.utils import flow_to_image
import matplotlib.pyplot as plt


def gen_of_array_raft(eui_array_file, xslice=None, yslice=None, device="cpu",
                      save_file=None):
    eui_array = np.load(eui_array_file)['eui_map_seq_crop_array']
    if xslice is None:
        xslice = slice(None)
    if yslice is None:
        yslice = slice(None)
    eui_array = eui_array[yslice, xslice,:].transpose(2,0,1)

    of_array_raft = np.zeros((eui_array.shape[0]-1, 2, eui_array.shape[1], eui_array.shape[2]), dtype=np.float64)

    eui_array = np.stack([eui_array, eui_array, eui_array], axis=1)
    eui_array = torch.from_numpy(eui_array)

    model = raft_large(pretrained=True).to(device)
    model = model.eval()

    

    for ii, (img1, img2) in enumerate(zip(eui_array[:,:,:], eui_array[1:,:,:])):
    # Note: it would be faster to predict batches of flows instead of individual flows
        img1 = preprocess(img1[None]).to(device)
        img2 = preprocess(img2[None]).to(device)

        list_of_flows = model(img1, img2)
        predicted_flow = list_of_flows[-1][0]
        of_array_raft[ii] = predicted_flow.detach().cpu().numpy()

        if ii == 180:
            flow_img = flow_to_image(predicted_flow)
            img = F.to_pil_image(flow_img.to("cpu"))
            plt.imshow(np.asarray(img))
            plt.show()


    if save_file is not None:
        np.savez(save_file, of_array_raft=of_array_raft)
    
def preprocess(batch):
    transforms = T.Compose(
        [
            T.ConvertImageDtype(torch.float32),
            T.Normalize(mean=0.5, std=0.5),  # map [0, 1] into [-1, 1]
        ]
    )
    batch = transforms(batch)
    return batch


if __name__ == '__main__':
    gen_of_array_raft('../../../sav/optical_flow/eui_map_crop_array_0500_0600_0670_0760_1024.npz', 
                      xslice=slice(0, 160), yslice=slice(0, 160), save_file='../../../sav/optical_flow/of_1024_east_1_raft.npz')




import torch, torch.nn as nn, torch.nn.functional as F, wandb
from torch.utils.data import DataLoader
from torchvision import datasets, transforms

dev = "cpu"

tf = transforms.ToTensor()
tr = DataLoader(datasets.MNIST("./data", train=True, download=True, transform=tf),
                batch_size=128, shuffle=True)
te = DataLoader(datasets.MNIST("./data", train=False, download=True, transform=tf),
                batch_size=512)

model = nn.Sequential(
    nn.Conv2d(1, 16, 3, padding=1), nn.ReLU(), nn.MaxPool2d(2),
    nn.Conv2d(16, 32, 3, padding=1), nn.ReLU(), nn.MaxPool2d(2),
    nn.Flatten(), nn.Linear(32 * 7 * 7, 10)).to(dev)

with wandb.init(project="cnn-demo",
                config={"epochs": 5, "lr": 1e-3, "batch_size": 128}) as run:
    opt = torch.optim.Adam(model.parameters(), lr=run.config.lr)
    run.define_metric("epoch")
    run.define_metric("test_accuracy", step_metric="epoch")
    run.watch(model, log="all", log_freq=100)

    for epoch in range(run.config.epochs):
        model.train()
        for xb, yb in tr:
            xb, yb = xb.to(dev), yb.to(dev)
            loss = F.cross_entropy(model(xb), yb)
            opt.zero_grad(); loss.backward(); opt.step()
            run.log({"train_loss": loss.item()})

        model.eval(); correct = 0
        with torch.no_grad():
            for xb, yb in te:
                correct += (model(xb.to(dev)).argmax(1).cpu() == yb).sum().item()

        run.log({"epoch": epoch + 1,
                 "test_accuracy": correct / len(te.dataset)})
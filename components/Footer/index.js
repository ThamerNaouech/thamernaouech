import React from "react";
import Socials from "../Socials";
import Link from "next/link";
import Button from "../Button";

const Footer = ({}) => {
  return (
    <>
      <div className="mt-5 laptop:mt-40 p-2 laptop:p-0">
        <div>
          <h1 className="text-2xl text-bold"></h1>
          <div className="mt-10">
            <h1 className="text-3xl tablet:text-2xl laptop:text-6xl laptopl:text-4xl text-bold">
              Let&apos;s get in
            </h1>
            <h1 className="text-3xl tablet:text-2xl laptop:text-6xl laptopl:text-4xl text-bold">
              touch
            </h1>
            <div className="mt-10">
              <Socials />
            </div>
          </div>
        </div>
      </div>
      <h1 className="text-sm text-bold mt-2 laptop:mt-10 p-2 laptop:p-0">
        Made With ❤ by{" "}
        <Link href="https://www.linkedin.com/in/thamernaouech/">
          <a>Thamer Naouech</a> 
        </Link>
      </h1>
    </>
  );
};

export default Footer;
